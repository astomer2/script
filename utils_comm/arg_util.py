import argparse
import logging
from datetime import datetime
from pathlib import Path

DATE_TIME = "%Y_%m_%d %H:%M:%S"


class ArgparseUtil(object):
    """
    参数解析工具类
    """

    def __init__(self):
        """Basic args"""
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("--seed", default=2, type=int)

    def af_batch_runner(self):
        """task args"""
        self.parser.add_argument(
            "--gpu_device_id", default=1, type=int, help="the GPU NO."
        )
        self.parser.add_argument(
            "--af_data", type=str, default="/mnt/sda/af/af_data", help=""
        )
        self.parser.add_argument(
            "--input_root_dir", type=str, default="/mnt/sda/af/af_input/tasks", help=""
        )
        self.parser.add_argument(
            "--out_root_dir", type=str, default="/mnt/sda/af/af_out/tasks", help=""
        )
        self.parser.add_argument(
            "--overwrite_predicted_pdb", type=int, default=0, help="0 false, 1 true"
        )
        self.parser.add_argument(
            "--read_full_cached_feature", type=int, default=1, help="0 false, 1 true"
        )
        self.parser.add_argument(
            "--model_preset", type=str, default="multimer", help=""
        )
        self.parser.add_argument("--task_name", type=str, default="", help="")
        self.parser.add_argument("--sub_task_name", type=str, default="", help="")
        self.parser.add_argument(
            "--run_relax", type=int, default=1, help="0 false, 1 true"
        )
        self.parser.add_argument(
            "--run_peptide_msa", type=int, default=1, help="0 false, 1 true"
        )
        self.parser.add_argument(
            "--force_peptide_re_feature", type=int, default=1, help="0 false, 1 true"
        )
        self.parser.add_argument(
            "--reverse", type=int, default=0, help="0 false, 1 true"
        )
        self.parser.add_argument(
            "--part_num", type=int, default=1, help="NB: starts from 1"
        )
        self.parser.add_argument(
            "--total_parts", type=int, default=2, help="total parts num of input files"
        )
        self.parser.add_argument(
            "--only_calc_seq_feature", type=int, default=0, help="0/1 bool"
        )
        self.parser.add_argument(
            "--only_run_model", type=int, default=0, help="0/1 bool"
        )
        self.parser.add_argument(
            "--collect_results",
            default=0,
            type=int,
            help="0 false, 1 true, not run alpha fold, but just collect results.",
        )
        self.parser.add_argument(
            "--merge_results",
            default=0,
            type=int,
            help="0 false, 1 true, not run alpha fold, but merge_results from different IPs.",
        )
        args = self.parser.parse_args()
        return args


def save_args(args, output_dir=".", with_time_at_filename=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    t0 = datetime.now().strftime(DATE_TIME)
    if with_time_at_filename:
        out_file = output_dir / f"args-{t0}.txt"
    else:
        out_file = output_dir / "args.txt"
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(f"{t0}\n")
        for arg, value in vars(args).items():
            f.write(f"{arg}: {value}\n")


def log_args(args, logger=None, save_dir=None):
    if logger is None:
        logger = logging.getLogger()
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(filename)s %(lineno)d: %(message)s",
            datefmt="%y-%m-%d %H:%M",
        )

    for arg, value in vars(args).items():
        logger.info(f"{arg}: {value}")

    if save_dir is not None:
        logger.info("Save args to the dir %s", save_dir)
        save_args(args, save_dir)
