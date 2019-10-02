echo running condenser
python condenser.py \
	samples/sample_sequences.fasta \
	--maxft 40 \
	--out sample_condenser_output.csv

echo rendering
python render.py \
	sample_condenser_output.csv \
	samples/sample_sequences.fasta \
	--region 130,469 \
	--alignment samples/HIV_env_alignment.csv \
	--phenotype samples/sample_phenotypes.csv \
	--out sample_renderer_output.pdf \
	--width 11 \
	-nonconserved
