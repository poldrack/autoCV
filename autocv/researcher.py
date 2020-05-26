"""
class for a researcher
"""

import os
import json


def get_params(params_file='params.json'):
    if os.path.exists(params_file):
        with open(params_file) as f:
            params = json.load(f)
    else:
        raise FileNotFoundError('Please create a json file called params.json containing the fields email (with your email address), orcid (with your ORCID id) and query (with your pubmed query)- see documentation for help')
    required_fields = ['address', 'lastname', 'firstname', 'email', 'orcid', 'query', 'url', 'phone']
    for field in required_fields:
        assert field in params
    return(params)

class Researcher:

    def __init__(self, params_file='params.json'):
        self.load_params(params_file)

    def load_params(self, params_file):
        if os.path.exists(params_file):
            with open(params_file) as f:
                params = json.load(f)
        else:
            raise FileNotFoundError("""Please create a json file called params.json 
                                       containing the fields email (with your email address), orcid (with your ORCID id) 
                                       and query (with your pubmed query)- see documentation for help')
                                       """)
        for field in params:
            setattr(self, field, params[field])

if __name__ == '__main__':
    r = Researcher('../tests/params.json')
    print(vars(r))
