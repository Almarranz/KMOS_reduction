# KMOS_reduction

##Rection fo Kmos data using esoReflex workflow.

The sky and the target are in different observation blocks,so we run esoReflex separetly in both block. 
There are two differents epochs, that esoReflex combine in a simgle block (COMBINE_obs). 

### WARNING!
I have added at some point on my life a fits file called RESPONSE_KKK.fits. It seems that you need this file for the pipeline
yield the results flux calibrated. I dont remember how I end up with this conclusio or where this file was located.

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
> When the standar star plot pops up set ***use_input_kernel = Fale** and **RERUN**

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

spectra_combines.py. Extracts and plot spectra from .fits files obtained 
diredctly with Qfitsview. 

cube_combine.py (is at iMac). It cut the different redcution according to the 
walength used to fit the spaectrum and put the in a single cube.