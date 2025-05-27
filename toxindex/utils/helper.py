# helper.py
import functools
import logging
import argparse

def handle_exceptions(catch_exceptions=(Exception,), default_value=None, log_errors=True):
    """
    Decorator factory to handle specified exceptions safely during function execution.

    It wraps a function in a try-except block, catches specified exceptions,
    optionally logs them, and returns a defined default value upon error.

    Args:
        catch_exceptions (tuple): A tuple of Exception types to catch.
                                  Defaults to (Exception,), catching all standard exceptions.
        default_value: The value to return if one of the specified exceptions is caught.
                       Defaults to None.
        log_errors (bool): If True, logs the caught exception details using the logging module.
                           Defaults to True. Requires logging to be configured.

    Returns:
        The decorator function.
    """
    def decorator(func):
        @functools.wraps(func) # Preserves function metadata (name, docstring, etc.)
        def wrapper(*args, **kwargs):
            try:
                # Execute the original function
                return func(*args, **kwargs)
            except catch_exceptions as e:
                # Handle the exception if caught
                if log_errors:
                    # Log the error with traceback information
                    logger.error(
                        f"Caught exception {type(e).__name__} in function '{func.__name__}': {e}",
                        exc_info=True # Set to False if traceback is not needed
                    )
                # Return the specified default value
                return default_value
            # Note: Exceptions not listed in catch_exceptions will propagate normally.
        return wrapper
    return decorator

# helper.py

import functools
import time
import os
import logging
import tempfile # To store lock file in a temporary directory

# --- Optional: Basic Logging Configuration ---
# Configure logging in your main script or use this default.
# logging.basicConfig(
#     level=logging.INFO,
#     format='%(asctime)s - %(levelname)s - %(name)s - %(message)s'
# )
# Get a logger instance specifically for this module
logger = logging.getLogger(__name__)

# --- Decorator Definition ---

def rate_limit_lockfile(seconds, lock_file_prefix="ratelimit_lock_"):
    """
    Decorator to rate-limit function execution based on a lock file timestamp.

    Prevents the decorated function from running if it has already run
    successfully within the specified `seconds` interval. Uses a file
    in the system's temporary directory to store the timestamp of the
    last successful execution.

    Args:
        seconds (int | float): The minimum time interval required between
                               successful executions.
        lock_file_prefix (str): A prefix for the lock file name. The function
                                name will be appended to this prefix.

    Returns:
        The decorator function.
    """
    def decorator(func):
        # Generate a unique lock file name based on the function
        # Use temp directory for better cross-platform compatibility and cleanup
        temp_dir = tempfile.gettempdir()
        lock_file_name = os.path.join(temp_dir, f"{lock_file_prefix}{func.__name__}.lock")

        @functools.wraps(func) # Preserves function metadata
        def wrapper(*args, **kwargs):
            current_time = time.time()
            last_run_time = 0.0

            # 1. Check the lock file for the last successful run time
            try:
                if os.path.exists(lock_file_name):
                    with open(lock_file_name, 'r') as f:
                        content = f.read().strip()
                        # Attempt to read the timestamp
                        if content:
                            last_run_time = float(content)
            except (IOError, ValueError) as e:
                # Handle issues reading the file or parsing the float
                logger.warning(
                    f"Could not read or parse timestamp from lock file '{lock_file_name}'. "
                    f"Proceeding as if interval has passed. Error: {e}"
                )
                last_run_time = 0.0 # Treat as if it hasn't run recently

            # 2. Decide whether to run the function
            elapsed_time = current_time - last_run_time
            if elapsed_time >= seconds:
                logger.debug(f"Interval of {seconds}s passed ({elapsed_time:.2f}s elapsed). Running '{func.__name__}'.")
                try:
                    # 3. Execute the original function
                    result = func(*args, **kwargs)

                    # 4. Update lock file ONLY on successful execution
                    try:
                        with open(lock_file_name, 'w') as f:
                            f.write(str(current_time))
                        logger.debug(f"Updated lock file '{lock_file_name}' with timestamp {current_time:.2f}.")
                    except IOError as e:
                        logger.error(f"Failed to update lock file '{lock_file_name}': {e}")

                    return result # Return the result of the original function

                except Exception as e:
                    # If the decorated function fails, log it but DO NOT update the lock file
                    logger.error(f"Function '{func.__name__}' raised an exception: {e}", exc_info=True)
                    # Re-raise the exception so the caller knows it failed
                    raise
            else:
                # Interval has not passed
                wait_time = seconds - elapsed_time
                logger.info(
                    f"Skipping '{func.__name__}'. Last run was {elapsed_time:.2f}s ago. "
                    f"Need to wait {wait_time:.2f}s more."
                )
                # Indicate that the function was skipped (e.g., return None or raise a specific exception)
                return None # Or: raise RateLimitExceededError("Rate limit exceeded")

        return wrapper
    return decorator

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')