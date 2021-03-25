import glob
import requests
import time
from datetime import datetime

while True: 
  paths = glob.glob('*.rda')
  if len(paths) >= 10: 
    req = requests.post('https://api.telegram.org/bot1019909271:AAGJNOFwMq0GfOCdoHmqH31rGB2OpC-mw9w/sendMessage', json={"text": "O2 10 RDA", "chat_id": "749625204"}, headers={'Content-type': 'application/json'})
  else: 
    print("T: " + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    time.sleep(120)
