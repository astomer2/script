import numpy as np
import torch
import random
import torch.optim as optim
import logging


logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO, datefmt='%y-%m-%d %H:%M',
    format='%(asctime)s %(filename)s %(lineno)d: %(message)s')


def set_seed(seed: int):
    """
    Helper function for reproducible behavior to set the seed in `random`, `numpy`, `torch` and/or `tf` (if installed).

    Args:
        seed (`int`): The seed to set.
    """
    logger.info('seed %s', seed)
    random.seed(seed)
    np.random.seed(seed)
    try:
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    except Exception as identifier:
        pass


def get_device(device_id=0):
    """ if device_id < 0: device = "cpu" """
    if device_id < 0:
        device = "cpu"
    elif torch.cuda.is_available():
        gpu_count = torch.cuda.device_count()
        if gpu_count > device_id:
            device = f"cuda:{device_id}"
        else:
            device = "cuda:0"
    else:
        device = "cpu"
    logger.info(f'device: {device}')
    return device


class CosineWarmupScheduler(optim.lr_scheduler._LRScheduler):
    def __init__(self, optimizer, warmup, max_iters):
        self.warmup = warmup
        self.max_num_iters = max_iters
        super().__init__(optimizer)

    def get_lr(self):
        lr_factor = self.get_lr_factor(epoch=self.last_epoch)
        return [base_lr * lr_factor for base_lr in self.base_lrs]

    def get_lr_factor(self, epoch):
        lr_factor = 0.5 * (1 + np.cos(np.pi * epoch / self.max_num_iters))
        if self.warmup > 0  and epoch <= self.warmup:
            lr_factor *= epoch * 1.0 / self.warmup
        return lr_factor


def nan_equal(a,b):
    try:
        np.testing.assert_equal(a,b)
    except AssertionError:
        return False
    return True


def models_are_equal(model1, model2):
    model1.vocabulary == model2.vocabulary
    model1.hidden_size == model2.hidden_size
    for a,b in zip(model1.model.parameters(), model2.model.parameters()):
        if nan_equal(a.detach().numpy(), b.detach().numpy()) == True:
            logger.info("true")


def run_compile_when_pt2(model, enable_complile=True, save_gpu_memory=True):
    """  """
    pt_version = torch.__version__.split('.')[0]
    if pt_version == '2' and enable_complile:
        if save_gpu_memory:
            _model = torch.compile(model)
        else:
            _model = torch.compile(model, mode='reduce-overhead')
    else:
        _model = model
    return _model


if __name__ == "__main__":
    get_device()