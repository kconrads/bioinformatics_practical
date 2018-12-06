import requests

GO_id = 'GO:0005634'
Aspect_url = 'http://ibeetle-base.uni-goettingen.de/ibp/go/GO2aspect/' + GO_id
Aspect_raw = requests.get(Aspect_url).text
Aspect_list = Aspect_raw.split(':')
Aspect_long = Aspect_list[2].replace('"', '').replace("}", '')
if Aspect_long == 'cellular_component':
    Aspect = 'C'
elif Aspect_long == 'biological_process':
    Aspect = 'P'
else: Aspect = 'F'

print(Aspect_raw)