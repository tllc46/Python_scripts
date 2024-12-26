from selenium import webdriver
from selenium.webdriver.common.by import By

driver=webdriver.Chrome()

driver.get(url="https://drive.google.com/drive/folders/1w_Fb9cx6LTfhClh4GroAtA_HXwLXmP0J")
input("after scroll down, press enter")

elem=driver.find_element(by=By.XPATH,value="/html/body/div[3]/div/div[3]/div[2]/div/div/c-wiz/div/c-wiz/div[1]/c-wiz/div[2]/c-wiz/div[1]/c-wiz/c-wiz/div")

elems=elem.find_elements(by=By.XPATH,value="*")

for i in elems:
    elem=i.find_element(by=By.XPATH,value="div")
    url=elem.get_attribute(name="data-id")
    url="https://drive.google.com/drive/folders/"+url
    elem=elem.find_element(by=By.XPATH,value="div[1]/div/div[2]/div/div")
    label=elem.get_attribute(name="aria-label")
    label=label.split(sep=":")[1].strip()
    print(label,url)

driver.close()
