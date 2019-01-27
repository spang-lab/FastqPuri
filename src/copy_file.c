#include <stdio.h>

int copy_file(char *old_filename, char  *new_filename)
{
  //fprintf(stderr, "Copy '%s' -> '%s.\n", old_filename, new_filename);
  FILE  *ptr_old, *ptr_new;
  int  a;

  if (!(ptr_old = fopen(old_filename, "rb"))) {
    fprintf(stderr, "Could not open file to be read \n");
    return -1;
  }

  if (!(ptr_new = fopen(new_filename, "wb"))) {
    fprintf(stderr, "Could not open file to be written to \n");
    return -1;
  }

  while(1)
  {
    a  =  fgetc(ptr_old);

    if(!feof(ptr_old))
      fputc(a, ptr_new);
    else
      break;
  }

  fclose(ptr_new);
  fclose(ptr_old);
  return  0;
}

