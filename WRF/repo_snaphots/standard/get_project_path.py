#!/bin/python3
from pathlib import Path
current_path=Path('.').resolve()
parts=current_path.parts
if parts[1]=='g':
    string=f'+gdata/{parts[3]}'
elif parts[1]=='scratch':
    string=f'+scratch/{parts[2]}'
else:
    string=''
print(string)
