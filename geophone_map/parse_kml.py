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
