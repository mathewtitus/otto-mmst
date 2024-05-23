FROM python:3

WORKDIR /Users/mtitus/Documents/GitHub/SUNSTRAND/otto-mmst/

COPY py/requirements.txt ./
RUN python3 -m pip install --no-cache-dir -r requirements.txt

COPY py/MMST.ipynb . 
CMD jupyter notebook --allow-root --ip 0.0.0.0 --no-browser . && netstat

