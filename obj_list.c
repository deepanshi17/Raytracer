#include "raytracer.h"

obj_list* ol_cons (object o, obj_list* os)
{
  if (os == NULL)
    return ol_singleton(o);
  obj_list* result = malloc(sizeof(obj_list));
  result->obj = o;
  result->next = os;
  return result;
}

obj_list* ol_singleton (object s)
{
  obj_list* single = malloc(sizeof(obj_list));
  single->obj = s;
  single->next = NULL;
  return single;
}

void ol_free (obj_list* ss)
{
  if (ss == NULL)
    return;
  else {
    ol_free(ss->next);
    free(ss);
  }
}

char *ol_tos(obj_list *os)
{
  /* FIXME protect against buffer overflow */
  char buf[2048];
  memset(buf,0,2048);
  strcat(buf,"(");
  while (os) {
    char *s = object_tos(os->obj);
    strcat(buf,s);
    if (os->next)
      strcat(buf,"\n");
    os=os->next;
  }
  strcat(buf,")");
  return strdup(buf);
}

void ol_show(FILE *f, obj_list *os)
{
  char *s = ol_tos(os);
  fprintf(f,"%s\n",s);
  free(s);
}

