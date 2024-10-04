from time import sleep

from selenium import webdriver
from selenium.webdriver.common.by import By

start_url="https://pbbs.kbs.co.kr/general/read.html?bbs_id=R2002-0281-03-240543&id=1789866&post_no=4206&page=1&post_header="
file=open(file="song_list_2.txt",mode="w",encoding="utf-8")
driver=webdriver.Chrome()
driver.get(url=start_url)

while True:
    try:
        sleep(1)
        elem=driver.find_element(by=By.XPATH,value="/html/body/div/div[2]/div/div[4]/div[1]")
        text=elem.text
        file.write(text+"\n")
        elem=driver.find_element(by=By.XPATH,value="/html/body/div/div[2]/div/div[4]/div[5]/ul/li[2]/div[3]/a")
        next_url=elem.get_attribute(name="href")
        driver.get(url=next_url)
    except:
        print("something went wrong")
        print(driver.current_url)
        file.close()
        break

driver.close()
file.close()
