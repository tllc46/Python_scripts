#####  input arguments  #####
#####   1. input file   #####
input_file="Gangneung.txt"
##### 2. fill index map #####
fill_map={"green":5,"blue":10,"yellow":3,"pink":7}
#############################

from time import sleep

import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

df=pd.read_csv(filepath_or_buffer=input_file,sep=" ",names=["lat","lon","fill","symbol","text"])

driver=webdriver.Chrome()
wait=WebDriverWait(driver=driver,timeout=20)

#네이버 로그인
driver.get(url="https://nid.naver.com/nidlogin.login")
input("after login, press enter")

#새 리스트
driver.get(url="https://pages.map.naver.com/save-pages/pc/all-list")
for key,value in fill_map.items():
    #"새 리스트 만들기" 버튼
    elem=driver.find_element(by=By.XPATH,value="/html/body/div/div/div/div[2]/button")
    elem.click()
    #"새 리스트명을 입력해주세요"
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[2]/div/div[2]/div[2]/div[1]/div/div/input")))
    elem.send_keys(key)
    #색상 선택 버튼
    elem=driver.find_element(by=By.XPATH,value=f"/html/body/div[2]/div/div[2]/div[2]/div[2]/div/button[{value}]")
    elem.click()
    #"완료" 버튼
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[2]/div/div[2]/div[3]/button")))
    elem.click()
    #리스트가 추가될 때까지 대기
    wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,f'/html/body/div/div/div/ul/li[button/div[2]/div[1]/span/text()="{key}"]')))

#네이버 지도
driver.get(url="https://map.naver.com/p")
for index,data in df.iterrows():
    latitude=data["lat"]
    longitude=data["lon"]
    fill=data["fill"]
    symbol=data["symbol"]
    text=data["text"]
    #"장소, 버스, 지하철, 도로 검색"
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div/div[1]/div/div/div/input")))
    elem.send_keys(f"{latitude} {longitude}",Keys.RETURN)
    #"저장 추가" 버튼
    try:
        elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div/div[2]/div[1]/div/div/div/div/div[1]/div[2]/div[2]/div[1]/button")))
    except:
        print(text,"| address doesn't exist, skip without saving")
        continue
    else:
        elem.click()
    #"메모, 별명, URL 추가" 버튼
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/div/button")))
    elem.click()
    #"지도 위에 표시될 별명을 남겨주세요"
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/div/ul/li[2]/div/div/input")))
    if symbol=="X":
        text+=" X"
    elem.send_keys(text)
    #리스트 선택
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,f'/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/ul/li[button/strong/text()="{fill}"]/button')))
    elem.click()
    #"저장" 버튼
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[3]/button")))
    elem.click()
    #장소가 추가될 때까지 대기
    wait.until(method=EC.text_to_be_present_in_element_attribute(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div[1]/div[2]/div[1]/div/div/div/div/div[1]/div[2]/div[2]/div[1]/button"),attribute_="aria-pressed",text_="true"))
    print(text,"| done saving")

driver.close()
