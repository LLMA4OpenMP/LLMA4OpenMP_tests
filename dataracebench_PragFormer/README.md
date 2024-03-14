# dataracebench by PragFormer

Since PragFormer does not generate code itself, there are no files here.

If you want to test PragFormer, go to

- Official:          [PragFormer Demo by Pragformer](https://huggingface.co/spaces/Pragformer/PragFormer-demo)
- Ours:              [PragFormer Demo by us](https://huggingface.co/spaces/xiaoming10086/PragFormer_demo)

Right now, the official demo has some problems and can't run. We found out the error message said that it cannot use the model on the internet (don't know why). So we downloaded all the models and put them in one folder, and it runs perfectly.

The results in this repo were tested using our demo, our demo will not be updated from now.

The commit our demo uses (mostly the newest except for PragFormer, we think there is a bug in the latest version, so we are using the second-newest. Time of writing 11-28-2023):

- [PragFormer](https://huggingface.co/Pragformer/PragFormer):                    	 9029422a43cb171afc4d6a7a4da8d1482bfe92b0

- [PragFormer_private](https://huggingface.co/Pragformer/PragFormer_private):       	f7d208992d735e0764507d438dcc5198deb7cfaf

- [PragFormer_reduction](https://huggingface.co/Pragformer/PragFormer_reduction):  	7ae5fb492b9f7cddc98d6da2fa244f4b7e18140d

Tip: use "run with docker" and run the demo on your own computer, it is much faster (but only if you are testing it heavily, the docker image is more than 5GB).
