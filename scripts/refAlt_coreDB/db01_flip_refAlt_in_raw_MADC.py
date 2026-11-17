#!/usr/bin/python3

def get_refAlt_flipped_markerIDs(flipped_markers):
	inp = open(flipped_markers, 'r')
	line = inp.readline()
	flipped = []
	while line:
		flipped.append(line.strip().split(',')[0])
		line = inp.readline()
	inp.close()
	return(flipped)


def flip_ref_alt(flipped, madc):
	inp = open(madc, 'r')
	outp = open(madc.replace('.csv', '_flipRefAlt.csv'), 'w')
	line = inp.readline()
	updated = []
	while line:
		line_array = line.strip().split(',')
		if line_array[1] in flipped:
			if 'Ref' in line_array[0]:
				line_array[0] = line_array[0].replace('Ref', 'Alt')
			elif 'Alt' in line_array[0]:
				line_array[0] = line_array[0].replace('Alt', 'Ref')
			outp.write(','.join(line_array) + '\n')
			updated.append(line_array[1])
		else:
			outp.write(line)
		line = inp.readline()
	print(f'\n# Updated {len(set(updated))} markers:')
	print(set(updated))
	inp.close()
	outp.close()


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('flipped_markers',
						help='Marker IDs with Ref and Alt alleles flipped')
	
	parser.add_argument('raw_madc',
						help='Raw MADC file to be updated')
	
	args = parser.parse_args()

	flipped = get_refAlt_flipped_markerIDs(args.flipped_markers)
	
	flip_ref_alt(flipped, args.raw_madc)
