def parse_kml():
    from pykml import parser

    namespace={"kml":"http://www.opengis.net/kml/2.2"}

    kml_file_name="Weonju.kml"
    kml_file=open(file=kml_file_name)
    doc=parser.parse(fileobject=kml_file)
    kml_file.close()
    root=doc.getroot()
    folder=root.Document.xpath("kml:Folder[kml:name='원주_지오폰_어레이_위치']",namespaces=namespace)[0]

    csv_file_name="Weonju.txt"
    csv_file=open(file=csv_file_name,mode="w")

    for placemark in folder.Placemark:
        name=placemark.name
        longitude=placemark.ExtendedData.xpath("kml:Data[@name='longitude(가안)']",namespaces=namespace)[0].value
        latitude=placemark.ExtendedData.xpath("kml:Data[@name='latitude(가안)']",namespaces=namespace)[0].value
        csv_file.write(f"{latitude} {longitude} {name}\n")

    csv_file.close()

def naver_map():
    import pandas as pd
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver.common.keys import Keys
    from selenium.webdriver.support.wait import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC

    list_name="강릉 출장"

    df=pd.read_csv(filepath_or_buffer="foo4.txt",sep=" ",names=["lat","lon","name","team"])

    driver=webdriver.Chrome()
    wait=WebDriverWait(driver=driver,timeout=10)

    #네이버 로그인
    driver.get(url="https://nid.naver.com/nidlogin.login")
    input("after login, press enter")

    #새 리스트
    driver.get(url="https://pages.map.naver.com/save-pages/pc/all-list")
    #"새 리스트 만들기" 버튼
    elem=driver.find_element(by=By.XPATH,value="/html/body/div/div/div/div[2]/button")
    elem.click()
    #"새 리스트명을 입력해주세요"
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[2]/div/div[2]/div[2]/div[1]/div/div/input")))
    elem.send_keys(list_name)
    #색상 선택 버튼
    color_id=12
    elem=driver.find_element(by=By.XPATH,value=f"/html/body/div[2]/div/div[2]/div[2]/div[2]/div/button[{color_id}]")
    elem.click()
    #"완료" 버튼
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[2]/div/div[2]/div[3]/button")))
    elem.click()
    wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,f'/html/body/div/div/div/ul/li[button/div[2]/div[1]/span/text()="{list_name}"]')))

    #네이버 지도
    driver.get(url="https://map.naver.com/p")
    for index,data in df.iterrows():
        latitude=data["lat"]
        longitude=data["lon"]
        station_name=data["name"]
        #"장소, 버스, 지하철, 도로 검색"
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[1]/div/div[2]/div[1]/div/div[1]/div/div/div/input")
        #elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div/div[1]/div/div/div/input")))
        elem.send_keys(f"{latitude} {longitude}")
        elem.send_keys(Keys.RETURN)
        #"저장 추가" 버튼
        try:
            elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div/div[2]/div[1]/div/div/div/div/div[1]/div[2]/div[2]/div[1]/button")))
        except:
            print(station_name,"| address doesn't exist, skip without saving")
            continue
        else:
            elem.click()
        #"메모, 별명, URL 추가" 버튼
        elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/div/button")))
        elem.click()
        #"지도 위에 표시될 별명을 남겨주세요"
        elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/div/ul/li[2]/div/div/input")))
        elem.send_keys(station_name)
        #리스트 선택
        elem=driver.find_element(by=By.XPATH,value=f'/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/ul/li[button/strong/text()="{list_name}"]/button')
        elem.click()
        #"저장" 버튼
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[3]/button")
        elem.click()
        wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[1]/div/div[2]/div[1]/div/div[2]/div[1]/div/div/div/div/div[1]/div[1]/div[2]/div[2]/div/button")))
        print(station_name,"| done saving")

    driver.close()

def kakao_map():
    import time

    import pandas as pd
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver.common.keys import Keys
    from selenium.webdriver.support.wait import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC

    list_name="강릉_출장"

    df=pd.read_csv(filepath_or_buffer="Gangneung.txt",sep=" ",names=["lat","lon","name","team"])
    color_map={"green":4,"blue":6,"yellow":2,"pink":7}

    driver=webdriver.Chrome()
    wait=WebDriverWait(driver=driver,timeout=10)

    #카카오 로그인
    driver.get(url="https://accounts.kakao.com/login/?continue=https%3A%2F%2Fmap.kakao.com")
    input("after login, press enter")

    #새 리스트
    #"지도 설정-교통저보, 지형도, 날씨 등 원하는 정보를 사용하세요."
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[10]")
    elem.click()
    #"MY" 버튼
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[5]/div[1]/div/div/ul/li[5]/a")
    elem.click()
    #"새 그룹 추가" 버튼
    elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[5]/div[5]/div/button[1]")))
    elem.click()
    #"그룹명을 입력하세요."
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[4]/form/fieldset/div[2]/dl[2]/dd/div/input")
    elem.send_keys(list_name)
    #"완료" 버튼
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[4]/form/fieldset/div[3]/button")
    elem.click()
    time.sleep(1)
    #"그룹이 생성되었습니다."
    driver.switch_to.alert.accept()

    for index,data in df.iterrows():
        latitude=data["lat"]
        longitude=data["lon"]
        station_name=data["name"]
        color_id=color_map[data["team"]]
        #"장소, 주소, 버스 검색"
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[2]/div/div/form/fieldset/div[1]/input")
        elem.clear()
        elem.send_keys(f"{latitude} {longitude}")
        elem.send_keys(Keys.RETURN)
        #즐겨찾기 추가
        elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[7]/div[6]/div[7]/div[2]/div/div[6]/div/div/div[2]/div[3]/div/div[1]/a[1]")))
        elem.click()
        #그룹 선택
        elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,f'/html/body/div[20]/div[2]/div[2]/ul/li[a/span[2]/strong/text()="{list_name}"]/a')))
        elem.click()
        #별명 입력
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[3]/form/fieldset/div[2]/div[1]/input")
        elem.click()
        elem.clear()
        elem.send_keys(station_name)
        #색상 선택 버튼
        elem=driver.find_element(by=By.XPATH,value=f"/html/body/div[20]/div[3]/form/fieldset/div[2]/ul/li[{color_id}]/input")
        elem.click()
        #"완료" 버튼
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[20]/div[3]/form/fieldset/div[3]/button")
        elem.click()
        #"즐겨찾기가 저장되었습니다."
        elem=wait.until(method=EC.visibility_of_element_located(locator=(By.XPATH,"/html/body/div[20]/div[7]/div[2]/a")))
        elem.click()
        print(station_name,"| done saving")
        time.sleep(1)

    driver.close()
