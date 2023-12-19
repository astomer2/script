import time
from DrissionPage import ChromiumPage
from DrissionPage import ChromiumOptions
from DrissionPage.common import By
import os
from tqdm import tqdm
# Create a custom download directory


def download_HPEP_result(download_directory, url_log):

    os.makedirs(download_directory, exist_ok=True)

    driver = ChromiumOptions().set_paths(browser_path="/snap/bin/chromium")
    driver.set_argument("--headless")
    page = ChromiumPage(addr_driver_opts=driver)
    page.download_set.by_DownloadKit()
    page.download_set.save_path(download_directory)

    with open(url_log) as f:
        url_list = f.readlines()
    url_list = [url.split()[0] for url in url_list if "%" not in url]
    
    # 使用tqdm 显示进度
    for url in tqdm(url_list):
        page.get(url)
        time.sleep(1)

        try:
            download_btn = page.ele((By.XPATH, "//BODY/CENTER[1]/TABLE[1]/TBODY[1]/TR[1]/TD[1]/DIV[1]/A[25]"))
            download_btn.click()
        except Exception as e:
            # Handle exception if download button not found or click failed
            print(f"Error downloading for URL {url}: {e}")
            continue
        time.sleep(5)
    page.close_tabs()
    page.quit()

if __name__ == "__main__":
    # 用argparse解析命令行
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-o","--download_directory", type=str, help="download directory")
    parser.add_argument("-i","--url_log", type=str, help="url_log path")
    args = parser.parse_args()
    download_directory = args.download_directory
    url_log = args.url_log

    download_HPEP_result(download_directory, url_log) 