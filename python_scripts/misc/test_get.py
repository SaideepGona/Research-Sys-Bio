import requests, json
import pprint

def get_request(accession):

    headers = {'accept': 'application/json'}
    url = "http://www.encodeproject.org/biosample/" + accession + "/?frame=object"
    response = requests.get(url, headers = headers)
    response_dict = response.json()

    return response_dict


def find_file_type(accession):
    response = get_request(accession)
    pp.pprint(response)
    file_type = response["file_type"]
    assembly = response["biological_replicates"]
    #merged_status = response[""]

    return (accession, file_type, assembly)

pp = pprint.PrettyPrinter(indent=4)
response = get_request("ENCSR000EWS")
pp.pprint(response)