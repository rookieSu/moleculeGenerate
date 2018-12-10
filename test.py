import re
str = "   12 1"

print(re.match(r"^\s+\d?",str).span()[1]-1)
