import requests

data = {
'subdata': 'MRIVHRLLKKLHSL' , 
'nmod': '' , 
'cmod': '' , 
'modifications[]': '' , 
'bridge': '' , 
'op': '提交' , 
}
cookie = 'GSAccessSession=1688356190325621; HASH_GSAccessSession=84E7E64DB15C1D0E91242EBDED0ACFF0674278E8; sensorsdata2015jssdkcross=%7B%22%24device_id%22%3A%2218919df53452d1-01868db54a749e-7e56547f-1600000-18919df5346bf5%22%7D; _ga_GN7NBTW165=deleted; Hm_lvt_b5bfec95c900ed8b81f8bbc5131e0f53=1688356213,1688537151,1689230777,1689645940; pageReferrInSession=https%3A//www.genscript.com.cn/peptide-handbook.html%3Fsrc%3Dleftbar; firstEnterUrlInSession=https%3A//www.genscript.com.cn/click_peptide_service.html; VisitorCapacity=1; PHPSESSID=dnjaoq4446b5m6182re9ct4q37; HASH_PHPSESSID=D5D4FDF57BB78D86A578C0F981DA00A3DCA30F26; sa_jssdk_2015_www=%7B%22distinct_id%22%3A%221688356190325621%22%2C%22props%22%3A%7B%22%24latest_referrer%22%3A%22%22%2C%22%24latest_referrer_host%22%3A%22%22%2C%22%24latest_traffic_source_type%22%3A%22%E7%9B%B4%E6%8E%A5%E6%B5%81%E9%87%8F%22%2C%22%24latest_search_keyword%22%3A%22%E6%9C%AA%E5%8F%96%E5%88%B0%E5%80%BC_%E7%9B%B4%E6%8E%A5%E6%89%93%E5%BC%80%22%7D%7D; _gid=GA1.3.1572938284.1689837901; SeqCustomer=1600423834%7C%E8%93%9D%E6%96%AF%E7%91%9E%7C2023-07-20%7CACTIVE; HASH_SeqCustomer=78A923A4402F7F5FFE184277D11F7FC899BD4D75; username=%E8%93%9D%E6%96%AF%E7%91%9E; HASH_username=CECEEAEFAC86FC62A422AF791B12765B19550EE0; GSID=1600423834; HASH_GSID=D6A0452F1C8C26534A131E9F72760DECD71F7EFA; live800_c_r=1689838163154_1689838253649_0_0_0; GSLastLogId=161492628; HASH_GSLastLogId=9ADDF8E51EEF083E05A542E86175912069F24488; Hm_lpvt_b5bfec95c900ed8b81f8bbc5131e0f53=1689843210; _gat=1; pt_s_1ba5fb37=vt=1689843210549&cad=; _ga=GA1.1.1869725139.1688356214; GSLogId=161492671; HASH_GSLogId=AF4ACFE44661936AC4B573087E0CD3C42A3730D6; pt_1ba5fb37=deviceId%3De1afab09-88e3-4274-a116-66030b90f91a%26sessionId%3D2bc7eaba-8f0f-4899-adad-4167aa270024%26accountId%3D%26vn%3D30%26pvn%3D1%26sact%3D1689843243366%26; _ga_GN7NBTW165=GS1.1.1689843210.8.0.1689843243.27.0.0'
headers = {'Cookie':cookie
}
a = requests.post(url='https://www.genscript.com.cn/tools/peptide-molecular-weight-calculator' , headers=headers , data=data)

with open('MRIVHRLLKKLHSL.html' , 'wb') as f:
    f.write(a.content)