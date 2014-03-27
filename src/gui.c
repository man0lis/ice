#include <stdio.h>
#include <stdlib.h>

#include "SDL.h"

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;

void putPixel(SDL_Renderer* renderer, int x, int y, int alpha) {
  SDL_SetRenderDrawColor(renderer, 0, 0, 0, alpha);
  SDL_RenderDrawPoint(renderer, x, y);
}

int main(int argc, char *argv) {

  SDL_Window* window = NULL;
  SDL_Renderer* renderer = NULL;
  SDL_Texture* texture = NULL;
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
      renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_TARGETTEXTURE);

      SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

      texture = SDL_CreateTexture(
          renderer,
          SDL_PIXELFORMAT_ARGB8888,
          SDL_TEXTUREACCESS_TARGET,
          SCREEN_WIDTH,
          SCREEN_HEIGHT
          );

      //clear screen
      SDL_SetRenderTarget(renderer, texture);
      SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
      SDL_RenderClear(renderer);
      SDL_SetRenderTarget(renderer, NULL);
      SDL_RenderCopy(renderer, texture, NULL, NULL);
      SDL_RenderPresent(renderer);

      while(!quit) {
        while(SDL_PollEvent(&event) != 0) {
          switch(event.type) {
            case SDL_QUIT:
                quit = 1;
                break;
            case SDL_MOUSEBUTTONUP:
                printf("Mouse up event\n");
                SDL_SetRenderTarget(renderer, texture);
                putPixel(renderer, event.motion.x, event.motion.y, 255);
                SDL_SetRenderTarget(renderer, NULL);
                SDL_RenderCopy(renderer, texture, NULL, NULL);
                SDL_RenderPresent(renderer);
                break;
            case SDL_MOUSEBUTTONDOWN:
                printf("Mouse down event\n");
                break;
            case SDL_MOUSEMOTION:
                printf("Mouse motion event\n");
                break;
          }
        }
      }

      SDL_DestroyWindow(window);
      SDL_Quit();


      return 0;

    }
  }

}
