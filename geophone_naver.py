import time

import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys

list_name="geophone test1"

df=pd.read_csv(filepath_or_buffer="/home/tllc46/Downloads/Milyang_test",sep=" ",names=["lat","lon","name"])

driver=webdriver.Chrome()

driver.get(url="https://nid.naver.com/nidlogin.login")
input("after login, press enter")

driver.get(url="https://pages.map.naver.com/save-pages/pc/all-list")
elem=driver.find_element(by=By.XPATH,value='//*[@id="app"]/div/div/div[2]/button')
elem.click()
time.sleep(secs=2)
elem=driver.find_element(by=By.XPATH,value='//*[@id="swt-save-input-folderview-list"]')
elem.send_keys(list_name)
color_id=10
elem=driver.find_element(by=By.XPATH,value=f'//*[@id="swt-save-widget-wrap"]/div[2]/div[2]/div[2]/div/button[{color_id}]')
elem.click()
elem=driver.find_element(by=By.XPATH,value='//*[@id="swt-save-widget-wrap"]/div[2]/div[3]/button')
elem.click()
time.sleep(secs=2)

driver.get(url="https://map.naver.com/p")
for i in range(len(df)):
    latitude=df.loc[i,"lat"]
    longitude=df.loc[i,"lon"]
    station_name=df.loc[i,"name"]
    elem=driver.find_element(by=By.CLASS_NAME,value="input_search")
    elem.send_keys(f"{latitude} {longitude}")
    elem.send_keys(Keys.RETURN)
    time.sleep(secs=2)
    elem=driver.find_element(by=By.XPATH,value='//*[@id="section_content"]/div/div/div/div/div[1]/div[2]/div[2]/div[1]/button')
    elem.click()
    time.sleep(secs=2)
    elem=driver.find_element(by=By.XPATH,value='//*[@id="swt-save-widget-wrap"]/div[2]/div[2]/div/button')
    elem.click()
    elem=driver.find_element(by=By.XPATH,value='//*[@id="swt-save-input-listview-nickname"]')
    elem.send_keys(station_name)
    elems=driver.find_elements(by=By.CLASS_NAME,value="swt-save-group-item")
    for elem in elems:
        sub_elem=elem.find_element(by=By.CLASS_NAME,value="swt-save-group-name")
        if sub_elem.text.split(sep="\n")[-1]==list_name:
            break
    elem=elem.find_element(by=By.CLASS_NAME,value="swt-save-group-info")
    elem.click()
    elem=driver.find_element(by=By.XPATH,value='//*[@id="swt-save-widget-wrap"]/div[2]/div[3]/button')
    elem.click()
    time.sleep(secs=2)

driver.close()
