import time
from DrissionPage import ChromiumPage
from DrissionPage import ChromiumOptions
from DrissionPage.common import By
import os 
# Create a custom download directory

download_directory = '/mnt/nas1/lanwei-125/PRLR/HPEP/test'
os.makedirs(download_directory, exist_ok=True)

url_log = "/mnt/nas1/lanwei-125/PRLR/HPEP/HPEP_dock/PRLR-global-log.txt"

if not os.path.exists(download_directory): 
    os.makedirs(download_directory)

driver = ChromiumOptions().set_paths(browser_path='/snap/bin/chromium')
driver.set_argument('--headless')
page = ChromiumPage(addr_driver_opts=driver)

page.download_set.by_DownloadKit()
page.download_set.save_path(download_directory)
# Read URLs from the log file

with open(url_log) as f:
    url_list = f.readlines()

# Filter out URLs containing time information
url_list = [url.split()[0] for url in url_list if '%' not in url]

for url in url_list:
    page.get(url)
    time.sleep(1)

    # Find the download button and click it
    download_btn = page.ele((By.XPATH, "//BODY/CENTER[1]/TABLE[1]/TBODY[1]/TR[1]/TD[1]/DIV[1]/A[25]"))
    #os.system("wget "+ download_btn)
    download_btn.click()

    # Wait for download to complete
    time.sleep(5)

# Close the current tab
page.close_tabs()

# Quit the browser
page.quit()
