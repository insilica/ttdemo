# simplecache.py

import os
import hashlib
import pickle
import functools

def simple_cache(cache_dir):
    os.makedirs(cache_dir, exist_ok=True)

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            key = f"{func.__name__}:{args}:{kwargs}"
            hash_key = hashlib.sha256(key.encode()).hexdigest()
            cache_path = os.path.join(cache_dir, f"{hash_key}.pkl")

            if os.path.exists(cache_path):
                with open(cache_path, 'rb') as f:
                    return pickle.load(f)

            result = func(*args, **kwargs)
            with open(cache_path, 'wb') as f:
                pickle.dump(result, f)
            return result

        return wrapper

    return decorator
