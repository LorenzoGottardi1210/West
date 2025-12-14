from setuptools import setup
import json

with open("../VERSION.json", "r") as file:
    data = json.load(file)

setup(version=data["version"])
