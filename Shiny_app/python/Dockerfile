# Limma DEG Explorer — Python Shiny (also used by the Hugging Face Space)
# Build and run from Shiny_app/python:
#   docker build -t limma-deg-explorer-py .
#   docker run --rm -p 7860:7860 limma-deg-explorer-py
#   open http://localhost:7860

FROM python:3.11-slim

# Hugging Face Spaces run containers as user 1000
RUN useradd -m -u 1000 user
WORKDIR /home/user/app

COPY --chown=user requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY --chown=user app.py analysis.py report.py ./
COPY --chown=user templates/ templates/
COPY --chown=user data/ data/

USER user
ENV HOME=/home/user MPLCONFIGDIR=/tmp/matplotlib
EXPOSE 7860
CMD ["shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "7860"]
