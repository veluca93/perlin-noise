from cperlin import PerlinNoiseFactory
import PIL.Image

size = 1200
res = 120
frames = 120
frameres = 40
space_range = size // res
frame_range = frames // frameres

pnf = PerlinNoiseFactory(
    3, octaves=4, tile=(space_range, space_range, frame_range), unbias=True)

for t in range(frames):
    img = PIL.Image.new('L', (size, size))
    # (scale, start, stop=start+1, step=1)
    data = pnf.get_range_as_bytes((res, 0, size), (res, 0, size),
                                  (frameres, t))
    img.putdata(data)

    img.save("noiseframe{:03d}.png".format(t))
    print(t)
