import time
from DrissionPage import ChromiumPage
from DrissionPage import ChromiumOptions
from DrissionPage.common import By
import os

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

    for url in url_list:
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
    download_directory = "/mnt/nas1/lanwei-125/CD44/docking/HPEP/"
    url_log = "/mnt/nas1/lanwei-125/CD44/docking/HPEP/OPN-log.txt"
    download_HPEP_result(download_directory, url_log) 