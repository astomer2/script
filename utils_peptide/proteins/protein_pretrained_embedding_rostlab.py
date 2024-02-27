import os
import sys
import re
import torch
from torch.nn.utils.rnn import pad_sequence
from transformers import T5EncoderModel, T5Tokenizer

sys.path.append(os.path.abspath('.'))

from icecream import ic

from utils_comm.log_util import logger
from utils_comm.train_util import get_device

ic.configureOutput(includeContext=True, argToStringFunction=str)

# esm model has limit on max seq len
MAX_PROTEIN_LEN = 2700


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

    def calc_embeddings_and_convert2cpu(self, pep_seqs, prot_seqs, categories):
        """  """
        with torch.no_grad():
            dataset = []
            pep_token_representations = self.__calc_embedding(pep_seqs)
            prot_token_representations = self.__calc_embedding(prot_seqs)
            for i, (pep_seq, prot_seq, category) in enumerate(zip(pep_seqs, prot_seqs, categories)):
                dataset.append([
                    pep_token_representations.last_hidden_state[i, :len(pep_seq)].cpu(),
                    prot_token_representations.last_hidden_state[i, :len(prot_seq)].cpu(),
                    category
                ])
        return dataset

    def __calc_embedding(self, seqs):
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

    def calc_embedding(self, seqs):
        """  """
        with torch.no_grad():
            token_representations = self.__calc_embedding(seqs)
            tensors = []
            for i, seq in enumerate(seqs):
                tensors.append(token_representations[i, : len(seq)])
            padded_tensor = pad_sequence(tensors, batch_first=True)
        return padded_tensor

    
def test_prot_transformer(input_seqs):
    """  """
    embedder = ProtTransEmbedding(gpu_id=0)
    logger.info('%s', embedder.tokenizer.get_vocab())
    result = embedder.calc_embedding(input_seqs)
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