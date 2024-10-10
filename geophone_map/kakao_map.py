#####   input arguments   #####
#####    1. input file    #####
input_file="Gangneung.txt"
#####  2. color index map #####
color_map={"green":4,"blue":6,"yellow":2,"pink":7}
##### 3. symbol index map #####
symbol_map={"O":1,"X":8}
###############################

from time import sleep

import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

df=pd.read_csv(filepath_or_buffer=input_file,sep=" ",names=["lat","lon","station","color","symbol"])

driver=webdriver.Chrome()
wait=WebDriverWait(driver=driver,timeout=10)
wait_x=WebDriverWait(driver=driver,timeout=2)

#카카오 로그인
driver.get(url="https://accounts.kakao.com/login/?continue=https%3A%2F%2Fmap.kakao.com")
input("after login, press enter")

#"지도 설정-교통저보, 지형도, 날씨 등 원하는 정보를 사용하세요."
elem=driver.find_element(by=By.XPATH,value="/html/body/div[10]")
elem.click()

#새 그룹
#"MY" 버튼
elem=driver.find_element(by=By.XPATH,value="/html/body/div[5]/div[1]/div/div/ul/li[5]/a")
elem.click()
for key,value in symbol_map.items():
    #"새 그룹 추가" 버튼
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[5]/div[5]/div/button[1]")))
    elem.click()
    #표식 설정
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[4]/form/fieldset/div[2]/div[1]/button")
    elem.click()
    #표식 결정
    elem=driver.find_element(by=By.XPATH,value=f"/html/body/div[20]/div[4]/form/fieldset/div[2]/div[1]/ul/li[{value}]/a")
    elem.click()
    #"그룹명을 입력하세요."
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[4]/form/fieldset/div[2]/dl[2]/dd/div/input")
    elem.send_keys(key)
    #"완료" 버튼
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[4]/form/fieldset/div[3]/button")
    elem.click()
    sleep(1)
    #"그룹이 생성되었습니다."
    driver.switch_to.alert.accept()

for index,data in df.iterrows():
    latitude=data["lat"]
    longitude=data["lon"]
    station_name=data["station"]
    color_idx=color_map[data["color"]]
    symbol_name=data["symbol"]
    #"장소, 주소, 버스 검색"
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[2]/div/div/form/fieldset/div[1]/input")
    elem.clear()
    elem.send_keys(f"{latitude} {longitude}",Keys.RETURN)
    #즐겨찾기 추가
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[7]/div[6]/div[7]/div[2]/div/div[6]/div/div/div[2]/div[3]/div/div[1]/a[1]")))
    elem.click()
    #그룹 선택
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,f'/html/body/div[20]/div[2]/div[2]/ul/li[a/span[2]/strong/text()="{symbol_name}"]/a')))
    elem.click()
    #별명 입력
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[3]/form/fieldset/div[2]/div[1]/input")
    elem.click()
    elem.clear()
    elem.send_keys(station_name)
    #색상 선택 버튼
    elem=driver.find_element(by=By.XPATH,value=f"/html/body/div[20]/div[3]/form/fieldset/div[2]/ul/li[{color_idx}]/input")
    elem.click()
    #"완료" 버튼
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[3]/form/fieldset/div[3]/button")
    elem.click()
    #"즐겨찾기가 저장되었습니다."
    try:
        elem=wait_x.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[20]/div[7]/div[2]/a")))
    except:
        pass
    else:
        elem.click()
    #wait.until(method=EC.text_to_be_present_in_element_attribute(locator=(By.XPATH,"/html/body/div[7]/div[6]/div[7]/div[2]/div/div[6]/div[2]/div/div[2]/div[3]/div/div[1]/a[1]"),attribute_="class",text_="fav ACTIVE"))
    print(station_name,"| done saving")

driver.close()
