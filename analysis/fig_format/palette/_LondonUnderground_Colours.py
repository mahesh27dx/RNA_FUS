#__all__ = ['build_palette']

color_dict = {
        'Bakerloo': '#894e24',
        'Central': '#dc241f',
        'Circle': '#ffce00',
        'District': '#007229',
        'Hammersmith': '#d799af',
        'Jubilee': '#868f98',
        'Metropolitan': '#751056',
        'Northern': '#000000',
        'Piccadilly': '#0010a8',
        'Victoria': '#00a0e2',
        'Waterloo': '#76d0bd',
        'Overground': '#ff6600',
        'Docklands': '#009999',
        'Tramlink': '#66cc00',
        'Elizabeth': '#614199',
    }

class build_palette:

    def __getitem__(self, key):
        return color_dict[key]

    # return colours in a defined order
    color_cycle = [color_dict[s] for s in [
        'Bakerloo',
        'Central',
        'Circle',
        'District',
        'Hammersmith',
        'Jubilee',
        'Metropolitan',
        'Northern',
        'Piccadilly',
        'Victoria',
        'Waterloo',
        'Overground',
        'Docklands',
        'Tramlink',
        'Elizabeth'
    ]]
