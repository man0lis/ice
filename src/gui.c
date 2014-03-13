#include <stdio.h>
#include <stdlib.h>

#include "SDL.h"

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;

int main(int argc, char *argv) {

  SDL_Window* window = NULL;
  SDL_Surface* screenSurface = NULL;
  SDL_Event event;

  int quit = 0;

  if(SDL_Init(SDL_INIT_VIDEO) < 0) {
    printf("SDL could not initialize! SDL_Error: %s\n", SDL_GetError());
  } else {
    window = SDL_CreateWindow(
        "SDL tutorial",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        SCREEN_WIDTH,
        SCREEN_HEIGHT,
        SDL_WINDOW_SHOWN
        );
    if(window == NULL) {
      printf("Window could not be created! SDL_Error: %s\n", SDL_GetError());
    } else {
      screenSurface = SDL_GetWindowSurface(window);
      SDL_FillRect(screenSurface, NULL, SDL_MapRGB(screenSurface->format, 0xFF, 0xFF, 0xFF));
      SDL_UpdateWindowSurface(window);

      while(!quit) {
        while(SDL_PollEvent(&event) != 0) {
          if(event.type == SDL_QUIT) {
            quit = 1;
          }
        }
      }

      SDL_DestroyWindow(window);
      SDL_Quit();


      return 0;

    }
  }

}
