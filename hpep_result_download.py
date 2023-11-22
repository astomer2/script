from selenium import webdriver
from selenium.webdriver.chrome.service import Service
import time
from selenium.webdriver.common.by import By
from webdriver_manager.chrome import ChromeDriverManager

'''
s = Service("/usr/bin/chromedriver")
options = webdriver.ChromeOptions()
options.add_argument('--headless')
options.add_argument('ignore-certificate-errors') 
options.binary_location = "/usr/bin/chromium-browser"
'''
driver_path = ChromeDriverManager().install()
service = Service(driver_path)
# 创建 Chrome WebDriver，并指定服务
browser = webdriver.Chrome(service=service)

#browser = webdriver.Chrome() # 打开Chrome浏览器

# 从txt文件读取url
url_log = r"C:\Users\123\Desktop\jobwork\subject\MC1R\screen\v3\hpep\global\MC1R-global-log.txt"
with open(url_log) as f:
    url_list = f.readlines()

# 过滤掉包含时间信息的url    
url_list = [url.split()[0] for url in url_list if '%' not in url]


for url in url_list:

    browser.get(url)
    
    # 等待页面加载完成
    time.sleep(1)  

    # 找到下载按钮,点击下载
    download_btn = browser.find_element(By.XPATH,"//BODY/CENTER[1]/TABLE[1]/TBODY[1]/TR[1]/TD[1]/DIV[1]/A[25]")
    download_btn.click()

    # 等待下载完成
    time.sleep(2.5)
    
# 关闭当前标签页
browser.close()

# 关闭浏览器
browser.quit()