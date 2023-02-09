# KMOS_reduction

##Rection fo Kmos data using esoReflex workflow.

The sky and the target are in different observation blocks,so we run esoReflex separetly in both block. 
There are two differents epochs, that esoReflex combine in a simgle block (COMBINE_obs). 
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

> Directly on the cambas (**Can´t remember why these ones**)
> 
>> **molefict/calctrans : false**
>> **telluric and response correction : 0**

2. Run for target block and sky block **separately**.
> You will get SCI_RECONSTRUCTED_object and SCI_RECONSTRUCTED_sky files.

3. Createa .sof file (this is a regular text file with extension .sof) in this way:
SCI_RECONSTRUCTED_object.fits OBJECT_CUBE
 
SCI_RECONSTRUCTED_sky.fits SKY_CUBE 
> NOTE: as far as I know you have to do this one by one. 
>Meaning the .sof only can have on object file and one sky file. 
>If it contains more they have to be commented.
4. Use **sky_tweak** algortim with the .sof file created in steep 3. **Remember**: only *one* OBJECT_CUBE and  *one* SKY_CUBE at the time
> NOTE: sky_teawk algortin give the same name each time, 
>be sure you change it each time you running it (somethig like Ob_sky_tweak_file1.fits and so on)
4. With **kmos_combine --method=header** combine the .fit files previously obtained in step 3 

Additional scripts to esoReflex KMOS workflow
