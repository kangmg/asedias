from datetime import datetime
import json
from io import StringIO
from typing import Union

def json_dump(resultDict:dict, job_name:str, metadata:Union[str, dict])->None:
    """
    Dump the DIAS results into a JSON file.
  
    Parameters
    ----------
    resultDict : dict
        The result dictionary containing the DIAS results.
    job_name : str
        The name of the job, which will be used as the filename.
    metadata : Union[str, dict]
        Metadata associated with the job, can be a string or a dictionary.
    """
    _savepath = f'./{job_name}.json'
    with open(_savepath, "w") as file:
        json.dump(
            {
                "METADATA": metadata,
                "RESULT": resultDict
            }, 
            fp=file, 
            indent=4, 
            ensure_ascii=False
            )