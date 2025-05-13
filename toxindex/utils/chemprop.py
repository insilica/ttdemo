import aiohttp
import asyncio
import requests
import time
from tenacity import retry, stop_after_attempt, wait_exponential
from functools import wraps

def make_safe(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            return {"result": None, "error": str(e)}
    return wrapper


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=2, max=10))
def get_chemprop_prediction(inchi: str, property_token: str) -> dict:
    base_url = "http://chemprop-transformer-alb-2126755060.us-east-1.elb.amazonaws.com/predict"
    params = {"property_token": property_token, "inchi": inchi}
    response = requests.get(base_url, params=params)
    response.raise_for_status()  # Raise exception for bad status codes
    return response.json()

def get_chemprop_prediction_safe(inchi: str, property_token: str, retries: int = 5, delay: int = 2) -> dict:
    return make_safe(get_chemprop_prediction)(inchi, property_token, retries, delay)

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def chemprop_predict_all(inchi: str) -> list[dict]:
    base_url = "http://chemprop-transformer-alb-2126755060.us-east-1.elb.amazonaws.com/predict_all"
    params = {"inchi": inchi}
    response = requests.get(base_url, params=params)
    response.raise_for_status()  # Raise exception for bad status codes
    return response.json()
