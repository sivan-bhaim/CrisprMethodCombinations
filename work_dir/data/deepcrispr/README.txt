Generating outputs for a DeepCRISPR comparison
==============================================
1. Install Docker
		sudo apt -y install docker.io
		sudo systemctl enable --now docker
		
2. Obtain the DeepCRISPR docker container:
		sudo docker pull michaelchuai/deepcrispr:latest
		
3. Run the container interactively, with the main work directory as a volume:
		sudo docker run --name deep -v $WORK_DIR:/work_dir --rm -i -t michaelchuai/deepcrispr bash

-- The following isntruction should be performed inside the container --

4. Copy the run_deepcrispr.py script to the DeepCRISPR directory:
		cp /work_dir/scripts/tool_utils/run_deepcrispr.py /root/DeepCRISPR/
		
5. For each of the desired datasets other than Mixture, run the sequence-only
   regression variant (where <dataset name> determines the output's name, so,
   for the two Xu datasets, <dataset dir>=xu while the names should be xu and
   nr_xu, depending on the chromosome):
		python /root/DeepCRISPR/run_deepcrispr.py -s /work_dir/data/<dataset dir>/chr<chromosome number>.txt -o /work_dir/ -d <dataset name> -r
   Repeat for the epigenetics regression variant:
		python /root/DeepCRISPR/run_deepcrispr.py -s /work_dir/data/<dataset dir>/chr<chromosome number>.txt -o /work_dir/ -d <dataset name> -r -e
	Repeat for the epigenetics classification variant:
		python /root/DeepCRISPR/run_deepcrispr.py -s /work_dir/data/<dataset dir>/chr<chromosome number>.txt -o /work_dir/ -d <dataset name> -e

6. Produce the same output for the Mixture dataset:
		python /root/DeepCRISPR/run_deepcrispr.py -t /work_dir/data/mixture/mixture_targets.txt -o /work_dir/ -d mixture -r
		python /root/DeepCRISPR/run_deepcrispr.py -t /work_dir/data/mixture/mixture_targets.txt -o /work_dir/ -d mixture -r -e
		python /root/DeepCRISPR/run_deepcrispr.py -t /work_dir/data/mixture/mixture_targets.txt -o /work_dir/ -d mixture -e
