FROM --platform=linux/amd64 continuumio/miniconda3

WORKDIR /app
# Create the environment:
COPY demOMEGA.yml /app/demOMEGA.yml
RUN conda env create -f /app/demOMEGA.yml

# Make the conda env active by default:
RUN echo "source activate demOMEGA" > ~/.bashrc
ENV PATH /opt/conda/envs/synoSCO/bin:$PATH

CMD ["/bin/sh", "-c", "bash"]