import os
import sys
from functools import cache
import re
import esm
import torch
from torch.nn.utils.rnn import pad_sequence
from transformers import T5EncoderModel, T5Tokenizer

sys.path.append(os.path.abspath('.'))

from icecream import ic

from utils_comm.log_util import logger
from utils_comm.train_util import get_device

ic.configureOutput(includeContext=True, argToStringFunction=str)
from utils_comm.file_util import FileUtil

# esm model has limit on max seq len
MAX_PROTEIN_LEN = 2700


class EsmEmbedding(object):
    """
    """
    def __init__(self, gpu_id=0, max_len=None, keep_cls_eos=False) -> None:
        self.device = get_device(gpu_id)
        self.batch_tokenizer, self.model, self.esm_layer_num = get_esm_model(
            max_len=max_len)
        self.max_len=max_len
        pt_version = torch.__version__.split('.')[0]
        if pt_version == '2':
            self.model = torch.compile(self.model)
        self.model.to(self.device)
        self.keep_cls_eos = keep_cls_eos
        if keep_cls_eos:
            self.start_i = 0
        else:
            self.start_i = 1

    def cal_embedding_pairs_and_convert2cpu(self, pep_seq_inputs, prot_seq_inputs):
        """
        Returns: a list of tensor which are on cpu
        """
        with torch.no_grad():
            dataset = []
            pep_token_repre = self.cal_embedding(pep_seq_inputs)
            prot_token_repre = self.cal_embedding(prot_seq_inputs)
            for i, (pep_seq, prot_seq) in enumerate(zip(pep_seq_inputs, prot_seq_inputs)):
                if self.max_len:
                    pep_repre = pep_token_repre.cpu()
                    prot_repre = prot_token_repre.cpu()
                else:
                    if self.keep_cls_eos:
                        pep_end_i = len(pep_seq) + 2
                        prot_end_i = len(prot_seq) + 2
                    else:
                        pep_end_i = len(pep_seq) + 1
                        prot_end_i = len(prot_seq) + 1
                    pep_repre = pep_token_repre[i, self.start_i: pep_end_i].cpu()
                    prot_repre = prot_token_repre[i, self.start_i: prot_end_i].cpu()
                dataset.append([pep_repre, prot_repre])
        return dataset

    def cal_embedding(self, pep_seq_inputs):
        """ keeps gpu if gpu is available """
        batch_tokens = self.batch_tokenizer(pep_seq_inputs)
        batch_tokens = batch_tokens.to(self.device)
        results = self.model(
            batch_tokens, repr_layers=[self.esm_layer_num], return_contacts=False)
        token_representations = results["representations"][self.esm_layer_num]
        return token_representations

    def cal_and_pad_embedding(self, seqs):
        """ After cal_embedding, the for each seq, each pad_idx aa has different 1d values, and so remove padded idx 
        embeddings at end of seq """
        with torch.no_grad():
            token_representations = self.cal_embedding(seqs)
            if self.max_len:
                return token_representations
            tensors = []
            for i, seq in enumerate(seqs):
                if self.keep_cls_eos:
                    prot_end_i = len(seq) + 2
                else:
                    prot_end_i = len(seq) + 1
                tensors.append(token_representations[i, self.start_i: prot_end_i])
            padded_tensor = pad_sequence(tensors, batch_first=True)
        return padded_tensor


class ESMBatchTokenizer(object):
    """Callable to convert an unprocessed (labels + strings) batch to a
    processed (labels + tensor) batch.
    """

    def __init__(self, alphabet, truncation_seq_length: int = None, max_len=50):
        self.alphabet = alphabet
        self.truncation_seq_length = truncation_seq_length
        self.max_len = max_len

    def __call__(self, seq_str_list: list[str]):
        batch_size = len(seq_str_list)
        seq_encoded_list = [self.alphabet.encode(seq_str) for seq_str in seq_str_list]
        if self.truncation_seq_length:
            seq_encoded_list = [seq_str[:self.truncation_seq_length] for seq_str in seq_encoded_list]
        if isinstance(self.max_len, int) and self.max_len > 0:
            max_len = self.max_len
        else:
            max_len = max(len(seq_encoded) for seq_encoded in seq_encoded_list)
        tokens = torch.empty(
            (
                batch_size,
                max_len + int(self.alphabet.prepend_bos) + int(self.alphabet.append_eos),
            ),
            dtype=torch.int64,
        )
        tokens.fill_(self.alphabet.padding_idx)

        for i, seq_encoded in enumerate(seq_encoded_list):
            if self.alphabet.prepend_bos:
                tokens[i, 0] = self.alphabet.cls_idx
            seq = torch.tensor(seq_encoded, dtype=torch.int64)
            tokens[
                i,
                int(self.alphabet.prepend_bos) : len(seq_encoded)
                + int(self.alphabet.prepend_bos),
            ] = seq
            if self.alphabet.append_eos:
                tokens[i, len(seq_encoded) + int(self.alphabet.prepend_bos)] = self.alphabet.eos_idx

        return tokens


@cache
def get_esm_model(model_level='650m', max_len=None):
    if model_level == '650m':
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        esm_layer_num = 33
    elif model_level == '3b':
        model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
        esm_layer_num = 36
    batch_tokenizer = ESMBatchTokenizer(alphabet, max_len=max_len)
    model.eval()
    return batch_tokenizer, model, esm_layer_num


class ProtTransEmbedding(object):
    """
    embedding length is 1024
    """
    def __init__(self, gpu_id=0) -> None:
        self.device = get_device(gpu_id)
        self.model_root = '/mnt/sda/models/Rostlab/prot_t5_xl_half_uniref50-enc'
        self.tokenizer = T5Tokenizer.from_pretrained(self.model_root, do_lower_case=False)
        # Load the model
        self.model = T5EncoderModel.from_pretrained(self.model_root).to(self.device)
        # only GPUs support half-precision currently; if you want to run on CPU use full-precision (not recommended, much slower)
        self.model.full() if self.device=='cpu' else self.model.half()
        self.model = self.model.eval()

    def cal_embeddings_and_convert2cpu(self, pep_seqs, prot_seqs, categories):
        """  """
        with torch.no_grad():
            dataset = []
            pep_token_representations = self.cal_embedding(pep_seqs)
            prot_token_representations = self.cal_embedding(prot_seqs)
            for i, (pep_seq, prot_seq, category) in enumerate(zip(pep_seqs, prot_seqs, categories)):
                dataset.append([
                    pep_token_representations.last_hidden_state[i, :len(pep_seq)].cpu(),
                    prot_token_representations.last_hidden_state[i, :len(prot_seq)].cpu(),
                    category
                ])
        return dataset

    def cal_embedding(self, seqs):
        """  """
        logger.info('%s', seqs)
        _seqs = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in seqs]
        ids = self.tokenizer.batch_encode_plus(_seqs, add_special_tokens=True, padding="longest")
        logger.info('ids %s', ids)
        input_ids = torch.tensor(ids['input_ids']).to(self.device)
        attention_mask = torch.tensor(ids['attention_mask']).to(self.device)
        result = self.model(input_ids=input_ids, attention_mask=attention_mask)
        embedding_repr = result.last_hidden_state
        embedding_repr = embedding_repr.to(torch.float32)
        return embedding_repr

    def cal_and_pad_embedding(self, seqs):
        """  """
        with torch.no_grad():
            token_representations = self.cal_embedding(seqs)
            tensors = []
            for i, seq in enumerate(seqs):
                tensors.append(token_representations[i, : len(seq)])
            padded_tensor = pad_sequence(tensors, batch_first=True)
        return padded_tensor
    

def test_esm(input_seqs):
    """  """
    embedder = EsmEmbedding(gpu_id=0)
    logger.info('%s', embedder.batch_tokenizer.alphabet)
    ic(embedder.batch_tokenizer.alphabet.prepend_bos)
    ic(embedder.batch_tokenizer.alphabet.append_eos)
    ic(embedder.batch_tokenizer.alphabet.padding_idx)
    ic(embedder.batch_tokenizer.alphabet.cls_idx)
    ic(embedder.batch_tokenizer.alphabet.eos_idx)
    result = embedder.cal_embedding(input_seqs)
    logger.info(f' {result.shape}')
    logger.info(f' {result[0]}')
    if result.shape[0] > 1:
        logger.info(f' {result[1]}')

    
def test_prot_transformer(input_seqs):
    """  """
    embedder = ProtTransEmbedding(gpu_id=0)
    logger.info('%s', embedder.tokenizer.get_vocab())
    result = embedder.cal_embedding(input_seqs)
    logger.info(f' {result.shape}')
    logger.info(f' {result[0]}')
    if result.shape[0] > 1:
        logger.info(f' {result[1]}')


if __name__ == "__main__":

    seqs = [
        'ADE',
        'RH',
        'ADEaAX',
    ]
    long_seqs = [''.join(['A'] * 2700)]
    # ic(len(long_seqs), long_seqs)
    # long_seqs = FileUtil.read_raw_text('1.txt')[:16]
    # lens = [len(seq) for seq in long_seqs]
    # ic(max(lens))

    test_prot_transformer(seqs)