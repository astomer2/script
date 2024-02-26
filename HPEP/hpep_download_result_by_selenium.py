      
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from webdriver_manager.chrome import ChromeDriverManager

import time
def run(url_log,download_dir):
    options = Options()
    prefs = {'profile.default_content_setting_popups':0, "download.default_directory": f"{download_dir}"}
    options.add_experimental_option("prefs", prefs)

    driver_path = ChromeDriverManager().install()
    service = Service(driver_path)
    browser = webdriver.Chrome(service=service, options=options)

    with open(url_log) as f:
        url_list = f.readlines() 
    url_list = [url.split()[0] for url in url_list if '%' not in url]

    for url in url_list:

        browser.get(url)
        time.sleep(1) 
        try:
            download_btn = browser.find_element(By.XPATH,"//BODY/CENTER[1]/TABLE[1]/TBODY[1]/TR[1]/TD[1]/DIV[1]/A[25]")
            download_btn.click()
        except Exception as e:
            print(f"Error downloading for URL {url}: {e}")
            continue
        time.sleep(2.5)
    browser.close()
    browser.quit()

if __name__ == '__main__':

    url_log = r"Y:\lanwei-125\MC5R\dock\HPEP\cpep_sequence\log.txt"
    download_dir= r"Y:\lanwei-125\MC5R\dock\HPEP\cpep_sequence"
    run(url_log,download_dir)