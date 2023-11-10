# KMOS_reduction

##Rection fo Kmos data using esoReflex workflow.

The sky and the target are in different observation blocks,so we run esoReflex separetly in both block. 
There are two differents epochs, that esoReflex combine in a simgle block (COMBINE_obs). 

### WARNING!
I have added at some point on my life a fits file called RESPONSE_KKK.fits. It seems that you need this file for the pipeline
yield the results flux calibrated. I dont remember how I ended up with this conclusio or where this file was located.

Afeter visiting CAB in Madrid we made some changes in the pipe line configuration:
> We ara going to use the standard star for the teluric correction
> We are going to selecte 3 different wavelenghts ranges in the adjustement of the standard star
>> this way the adjusmente will be bettwer in the different zones of the spectrum
___
>NOTE: 
>Each OBs block is divided in two halfs, in between wich the sky was taken. That
>mean that the telescope moves and then alignment problem could arise when 
>combined the two halfs. So, it woudl be convinient to reduce only one of the 
>halfs and see the difference between reducing one and the two of then.
__

1. Change the following parametres:
> **Cube_reconstruction** actor on the canvas:
>
>> **init_no_substract = TRUE**
>>
>> **init_sky_tweak = FALSE

> Directly on the cambas (**Uses standard staras for teluric correction**)
> 
>> **molefict/calctrans : TRUE**
>> **telluric and response correction : 2**

#### OPTIONAL changes in the canvas:
The following changes could  improve the quality of both, the image and the spectrum

> number of wavelenght samples = 3072 (2048 por defecto)
> gobal_pixel_scale = 0.1 (0.2 por defecto) 

During the run of the pipeline:
> When the standar star plot pops up set ***use_input_kernel = False** and **RERUN**

2. Run for target block and sky block **separately**.
> You will get SCI_RECONSTRUCTED_object and SCI_RECONSTRUCTED_sky files.

3. sky_tweak.py. Ondce you have kmos run over science and sky cubes, this scripts
create the sof files with the SCI_RECONSTRUCTED and SKY_RECONTRUCTED files, 
runs the  **sky_tweak** algortim. Then crates a sof file with cubes subtracted 
by the sky and the runs **kmos_combine --method=header**  on it to generate 
the combined image.
> Note: you have to specify the folders with the scicience and sky cubes are 
> in the script. Also, specify if you have **singles_cubes** or ** reconstructed_cubes**

# Other scripts
## SPECTRA

First we combine epcoh p105 and p107 in order to incrase the SNR. Since the
images from both epoch are not well alingned, we are goin to aling IFU by IFU
using just one slice from the cube. Then we apply the alingment to the whole
cube and calculate the mean cube. 
Then, manually we extract the spectra using **Qfitsview**

1. kmos_chop.py. Separates the ifu from the rest of the image. Then using ds9
you have to aling then. The scrip asks you to record the offsets of the 
x and y coordinates. When you are satisfied with the alingement run it again 
and save them.

2. ifu_combine.py. Uses the offsets parametres previously extracted for the 
alingment of the whole cube (for a single IFU)

3. spect_extract.py Extract the sky and the extracts the sky subtracted spectra
Firsr, choose the *Sky* option in the console, then you have to selecte 
5 areas free of stars. Then, you run the script again 
selecting the *Sci* option and click on the stars. A fits file will be saved for 
each spectrum
4. spec_comp.py. Plot the spectra and generates a region file with the id of 
each star and its SNR.

