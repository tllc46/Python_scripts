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
    import time

    import pandas as pd
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver.common.keys import Keys

    list_name="원주 출장"

    df=pd.read_csv(filepath_or_buffer="/home/tllc46/Downloads/Weonju.txt",sep=" ",names=["lat","lon","name"])

    driver=webdriver.Chrome()

    #네이버 로그인
    driver.get(url="https://nid.naver.com/nidlogin.login")
    input("after login, press enter")

    #새 리스트
    driver.get(url="https://pages.map.naver.com/save-pages/pc/all-list")
    #"새 리스트 만들기" 버튼
    elem=driver.find_element(by=By.XPATH,value="/html/body/div/div/div/div[2]/button")
    elem.click()
    time.sleep(2)
    #"새 리스트명을 입력해주세요"
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[2]/div/div[2]/div[2]/div[1]/div/div/input")
    elem.send_keys(list_name)
    #색상 선택 버튼
    color_id=10
    elem=driver.find_element(by=By.XPATH,value=f"/html/body/div[2]/div/div[2]/div[2]/div[2]/div/button[{color_id}]")
    elem.click()
    #"완료" 버튼
    elem=driver.find_element(by=By.XPATH,value="/html/body/div[2]/div/div[2]/div[3]/button")
    elem.click()
    time.sleep(2)

    #네이버 지도
    driver.get(url="https://map.naver.com/p")
    for i in range(len(df)):
        latitude=df.loc[i,"lat"]
        longitude=df.loc[i,"lon"]
        station_name=df.loc[i,"name"]
        #"장소, 버스, 지하철, 도로 검색"
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[1]/div/div[2]/div[1]/div/div[1]/div/div/div/input")
        elem.send_keys(f"{latitude} {longitude}")
        elem.send_keys(Keys.RETURN)
        time.sleep(2)
        #"저장 추가" 버튼
        try:
            elem=driver.find_element(by=By.XPATH,value="/html/body/div[1]/div/div[2]/div[1]/div/div[2]/div[1]/div/div/div/div/div[1]/div[2]/div[2]/div[1]/button")
        except:
            print(station_name,"| address doesn't exist, skip without saving")
            continue
        else:
            elem.click()
            time.sleep(2)
        #"메모, 별명, URL 추가" 버튼
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/div/button")
        elem.click()
        #"지도 위에 표시될 별명을 남겨주세요"
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/div/ul/li[2]/div/div/input")
        elem.send_keys(station_name)
        #리스트 선택
        elem=driver.find_element(by=By.XPATH,value=f'/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[2]/ul/li[button/strong/text()="{list_name}"]/button')
        elem.click()
        #"저장" 버튼
        elem=driver.find_element(by=By.XPATH,value="/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div/div[2]/div[3]/button")
        elem.click()
        print(station_name,"| done saving")
        time.sleep(2)

    driver.close()

parse_kml()
naver_map()
