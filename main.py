"""Entry point notice.

This is a Streamlit application. To run it:

    uv run streamlit run app.py

Or with plain Python (after installing dependencies):

    streamlit run app.py

See README.md for full setup instructions.
"""

if __name__ == "__main__":
    import subprocess
    import sys
    subprocess.run([sys.executable, "-m", "streamlit", "run", "app.py"])
