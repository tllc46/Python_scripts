import pandas as pd
from lxml import etree

def parse_kml():
    kml_file="buan.kml"
    csv_file="buan.csv"

    namespace={None:"http://www.opengis.net/kml/2.2"}

    file=open(file=csv_file,mode="w")

    tree=etree.parse(source=kml_file)
    kml=tree.getroot()
    document=kml.find(path="Document",namespaces=namespace)

    for placemark in document.iterfind(path="Placemark",namespaces=namespace):
        name=placemark.find(path="name",namespaces=namespace).text
        description=placemark.find(path="description",namespaces=namespace)
        if description is not None:
            description=description.text
        else:
            description="no description"
        styleurl=placemark.find(path="styleUrl",namespaces=namespace).text[:17]
        point=placemark.find(path="Point",namespaces=namespace)
        coordinates=point.find(path="coordinates",namespaces=namespace).text
        coordinates=coordinates.strip().split(sep=",")
        print(coordinates[0],coordinates[1],coordinates[2],name,styleurl,description,sep=";",file=file)

    file.close()

def gen_kml():
    csv_file="buan.csv"
    kml_file="buan_restore.kml"

    style_list=["icon-1511-9C27B0","icon-1626-757575","icon-1899-000000","icon-1899-006064","icon-1899-0288D1","icon-1899-795548","icon-1899-7CB342","icon-1899-9C27B0","icon-1899-BDBDBD","icon-1899-C2185B","icon-1899-F57C00","icon-1899-FF5252","icon-1899-FFEA00"]
    color_list=["ffb0279c","ff757575","ff000000","ff646000","ffd18802","ff485579","ff42b37c","ffb0279c","ffbdbdbd","ff5b18c2","ff007cf5","ff5252ff","ff00eaff"]
    style_map=dict(zip(style_list,color_list))

    nsmap={None:"http://www.opengis.net/kml/2.2"}

    df=pd.read_csv(filepath_or_buffer=csv_file,sep=";",names=["lon","lat","elv","name","style","desc"])

    kml=etree.Element("kml",nsmap=nsmap)
    document=etree.SubElement(kml,"Document")
    etree.SubElement(document,"name").text="restored_map"

    for i in style_map:
        color_text=style_map[i]
        style=etree.SubElement(document,"Style",id=i)
        iconstyle=etree.SubElement(style,"IconStyle")
        etree.SubElement(iconstyle,"color").text=color_text
        icon=etree.SubElement(iconstyle,"Icon")
        etree.SubElement(icon,"href").text="https://www.gstatic.com/mapspro/images/stock/503-wht-blank_maps.png"

    for i,row in df.iterrows():
        placemark=etree.SubElement(document,"Placemark")
        etree.SubElement(placemark,"name").text=row["name"]
        etree.SubElement(placemark,"styleUrl").text=row["style"]
        extendeddata=etree.SubElement(placemark,"ExtendedData")
        data=etree.SubElement(extendeddata,"Data",name="위도")
        etree.SubElement(data,"value").text=str(row["lat"])
        data=etree.SubElement(extendeddata,"Data",name="경도")
        etree.SubElement(data,"value").text=str(row["lon"])
        data=etree.SubElement(extendeddata,"Data",name="설명")
        etree.SubElement(data,"value").text=row["desc"]
        point=etree.SubElement(placemark,"Point")
        etree.SubElement(point,"coordinates").text=str(row["lon"])+","+str(row["lat"])+","+str(row["elv"])

    tree=etree.ElementTree(element=kml)
    tree.write(file=kml_file,encoding="UTF-8",pretty_print=True,xml_declaration=True)
