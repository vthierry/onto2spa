# Here are how we compile the report

# Here are also a few links:
# - The publication repository page https://hal.inria.fr/hal-03338721
# - Access to the source of the code :
#   - The synthetical data : https://colab.research.google.com/drive/1Dulb2iwzDjqgS5WhTH_e-DCXAH_cprxP?usp=sharing
#     - Some of the code is changed and improved in the second notebook, not sure if that's very interesting
#   - Working with MNIST : https://colab.research.google.com/drive/1wklrZI3c5niIxBW3xgqQiB8AdeZ5nR17?usp=sharing
#     - There is a place where two functions are called, create_pi_model () and create_pi_conv_model (). The first gives a functional model, the second does not. 
#


all :
	rm -f *.log *.aux *.out *.toc 
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main
	rm *.log *.aux *.out *.toc 
