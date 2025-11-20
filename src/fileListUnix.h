/***** fileListUnix.h ******************************************************
 * Description: get file names from a directory 
 * Author: Mirjana Domazet-Loso, September 22, 2008
 *
 *****************************************************************************/ 
#ifndef FILELISTUNIX_H
#define FILELISTUNIX_H

struct dirent;
extern        int scandir(const char *restrict dirp,
			  struct dirent ***restrict namelist,
			  int (*filter)(const struct dirent *),
			  int (*compar)(const struct dirent **,
					const struct dirent **));

extern  int alphasort();

struct direct;
int selectFile(const struct direct *entry); /* pointer to this function is needed for listFilesUnix */
char **listFilesUnix(char *dir, int *count);

#endif /* FILELISTUNIX_H */

