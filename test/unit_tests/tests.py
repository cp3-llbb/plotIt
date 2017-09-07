from __future__ import division

import os
import re
import unittest
import shutil
import yaml
import tempfile
import subprocess

from configuration import get_configuration

class TemporaryFolder:
    def __init__(self):
        self.name = tempfile.mkdtemp()

    def __del__(self):
        shutil.rmtree(self.name, ignore_errors=True)

def get_golden_file(f):
    return os.path.join('golden', f)

convert_line_regexp = re.compile('(\d+):\s+\(\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)\)')
def get_images_likelihood(image1, image2):
    import subprocess

    compare = subprocess.Popen(['compare', image1, image2, '-compose', 'src', 'miff:-'], stdout=subprocess.PIPE)
    convert = subprocess.Popen(['convert', '-', '-depth', '8', '-define', 'histogram:unique-colors=true', '-format', '%[width] %[height]\n%c', 'histogram:info:-'], stdin=compare.stdout, stdout=subprocess.PIPE)

    compare.stdout.close()

    output = convert.communicate()[0]

    lines = output.split('\n')

    width, height = lines[0].split()

    pixels = int(width) * int(height)
    pixel_count = 0

    non_black_pixels = 0

    for i in xrange(1, len(lines)):
        matches = convert_line_regexp.search(lines[i])
        if matches is None:
            continue

        line_pixel = int(matches.group(1))
        pixel_count += line_pixel

        r, g, b, a = matches.group(2), matches.group(3), matches.group(4), matches.group(5)

        if int(r) == 255 and int(g) == 255 and int(b) == 255:
            # Black, ignore
            continue

        non_black_pixels += line_pixel

        if pixel_count == pixels:
            break

    return 1 - non_black_pixels / pixels

def to_png(image):
    """
    Convert image to PNG using 'convert'
    """
    import subprocess

    temp_name = None
    with tempfile.NamedTemporaryFile(delete=False) as f:
        temp_name = f.name + '.png'

    args = ['convert', '-define', 'pdf:use-cropbox=true', '-units', 'PixelsPerInch', '-density', '120', image, temp_name]
    subprocess.check_call(args)

    return temp_name


def upload_images(image1, image2):
    import subprocess

    diff_output = None
    with tempfile.NamedTemporaryFile(delete=False) as f:
        diff_output = f.name + '.pdf'

    # Create nice diff of image1 and image2
    subprocess.call(['compare', '-define', 'pdf:use-cropbox=true', image1, image2, diff_output])

    # Convert images to png
    image1 = to_png(image1)
    image2 = to_png(image2)
    diff_output = to_png(diff_output)

    # Upload to imgur
    from imgurpython import ImgurClient

    client_id = 'e32a7007635bb1f'
    client = ImgurClient(client_id, None)

    uploaded_image1 = client.upload_from_path(image1)
    uploaded_image2 = client.upload_from_path(image2)
    uploaded_diff = client.upload_from_path(diff_output)

    print "\nImages uploaded:"
    print " - Reference image: %s" % uploaded_image2['link']
    print " - Generated image: %s" % uploaded_image1['link']
    print " - Difference: %s" % uploaded_diff['link']

    print "\n-----"
    print "To delete image, execute these commands:"

    curl_command = 'curl -X DELETE -H "Authorization: Client-ID %s" https://api.imgur.com/3/image/%%s' % (client_id)

    print curl_command % uploaded_image1['deletehash']
    print curl_command % uploaded_image2['deletehash']
    print curl_command % uploaded_diff['deletehash']
    print "-----"
    print ""

class plotItSimpleTestCase(unittest.TestCase):
    def __init__(self, methodName='runTest'):
        super(plotItSimpleTestCase, self).__init__(methodName)

        # Switch to True to generate golden images
        self.__generate_golden_images = False

    def run_plotit(self, configuration):
        with tempfile.NamedTemporaryFile() as yml:
            yml.write(yaml.dump(configuration))
            yml.flush()
            with open(os.devnull, 'w+b') as null:
                subprocess.check_call(['../plotIt', yml.name, '-o', self.output_folder.name], stdout=null)

    def setUp(self):
        self.output_folder = TemporaryFolder()

    def tearDown(self):
        del self.output_folder

    def compare_images(self, image1, image2, threshold=0.995):
        if self.__generate_golden_images:
            # Just copy image1 to image2
            shutil.copyfile(image1, image2)
            print("Golden image %s generated" % image2)
            return

        likelihood = get_images_likelihood(image1, image2)

        try:
            self.assertGreater(likelihood, threshold, "Images too different: %.2f%% of similitude" % (likelihood * 100))
        except AssertionError as e:
            # Upload both image1, image2 and the diff
            # and output urls
            upload_images(image1, image2)

            raise e

    def suite(self):
        return unittest.TestLoader().loadTestsFromTestCase(self.__class__)


class plotItTestCase(plotItSimpleTestCase):
    def test_default_no_ratio(self):
        configuration = get_configuration()

        configuration['plots']['histo1']['show-ratio'] = False

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_no_ratio.pdf')
                )

    def test_default_ratio(self):
        configuration = get_configuration()

        configuration['plots']['histo1']['show-ratio'] = True

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_ratio.pdf')
                )

    def test_legend(self):
        configuration = get_configuration()

        configuration['legend']['columns'] = 1

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_1column_legend.pdf')
                )

        configuration['legend']['columns'] = 2

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_2columns_legend.pdf')
                )

        configuration['legend']['columns'] = 3

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_3columns_legend.pdf')
                )

        configuration['legend']['columns'] = 2
        configuration['files']['MC_sample1.root']['legend-order'] = 1
        configuration['files']['MC_sample2.root']['legend-order'] = 0

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_2columns_samplesordering_legend.pdf')
                )

        
        configuration = get_configuration()
        configuration['legend']['columns'] = 2
        configuration['legend']['position'][2] = 0.75

        configuration['legend']['entries'] = [
                {'label': 'Extra MC legend entry', 'type': 'mc', 'order': 1000},
                {'label': 'Extra signal legend entry', 'type': 'signal', 'order': 1000},
                ]


        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_extra_legend_entries.pdf')
                )

    def test_lines(self):
        configuration = get_configuration()

        configuration['configuration']['line-color'] = "#1685f2"
        configuration['configuration']['line-width'] = 3
        configuration['configuration']['line-type'] = 4

        configuration['plots']['histo1']['vertical-lines'] = [3, 7]
        configuration['plots']['histo1']['horizontal-lines'] = [80]
        configuration['plots']['histo1']['lines'] = [
                [[3, 20], [8, 90]],
                [[5.5, float('nan')], [5.5, float('nan')]]
                ]

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_lines.pdf')
                )

        configuration = get_configuration()

        configuration['configuration']['line-color'] = "#1685f2"
        configuration['configuration']['line-width'] = 3
        configuration['configuration']['line-type'] = 4

        configuration['plots']['histo1']['vertical-lines'] = [{
                'line-color': '#ae4fea',
                'line-width': 10,
                'value': 3}, 7]
        configuration['plots']['histo1']['horizontal-lines'] = [80]
        configuration['plots']['histo1']['lines'] = [
                [[3, 20], [8, 90]],
                [[5.5, float('nan')], [5.5, float('nan')]]
                ]

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_lines_inline_style.pdf')
                )

    def test_systematics(self):
        configuration = get_configuration()

        configuration['configuration']['luminosity-error'] = 0.4

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_lumi_error_0p4.pdf')
                )

        configuration = get_configuration()

        configuration['configuration']['luminosity-error'] = 0.
        configuration['systematics'] = []
        configuration['systematics'] += [{'syst1': 1.4}]

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_one_syst_0p4.pdf')
                )

        configuration = get_configuration()

        configuration['configuration']['luminosity-error'] = 0.
        configuration['systematics'] = []
        configuration['systematics'] += [{'syst1': 1.2}]
        configuration['systematics'] += [{'syst2': {'type': 'ln', 'prior': 1.2}}]

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_two_systs_0p28.pdf')
                )

        configuration = get_configuration()

        configuration['configuration']['luminosity-error'] = 0.
        configuration['systematics'] = ['alpha', 'beta']

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_two_systs_shape.pdf')
                )

        configuration = get_configuration()

        configuration['configuration']['luminosity-error'] = 0.
        configuration['systematics'] = ['i-do-not-exist']

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_syst_not_found.pdf')
                )

    def test_blinded(self):
        configuration = get_configuration()
        configuration['plots']['histo1']['blinded-range'] = [3, 5.2]

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_blinded_range.pdf')
                )

        configuration = get_configuration()

        configuration['plots']['histo1']['blinded-range'] = [3, 5.2]
        configuration['plots']['histo1']['log-y'] = True

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1_logy.pdf'),
                get_golden_file('default_configuration_blinded_range_logy.pdf')
                )

    def test_multi_stacks(self):
        configuration = get_configuration()

        # Remove data
        del configuration['files']['data.root']

        configuration['files']['MC_sample1.root']['stack-index'] = 0
        configuration['files']['MC_sample2.root']['stack-index'] = 1

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_multi_stacks.pdf')
                )

        configuration['files']['MC_sample1.root']['stack-index'] = 1
        configuration['files']['MC_sample2.root']['stack-index'] = 0

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_multi_stacks_reverse_ordering.pdf')
                )

        configuration['files']['MC_sample1.root']['stack-index'] = 0
        configuration['files']['MC_sample2.root']['stack-index'] = 1

        configuration['files']['MC_sample1.root']['fill-type'] = 0
        configuration['files']['MC_sample1.root']['line-type'] = 1
        configuration['files']['MC_sample1.root']['line-width'] = 2
        configuration['files']['MC_sample1.root']['legend-style'] = 'l'

        configuration['files']['MC_sample2.root']['fill-type'] = 0
        configuration['files']['MC_sample2.root']['line-type'] = 1
        configuration['files']['MC_sample2.root']['line-width'] = 2
        configuration['files']['MC_sample2.root']['legend-style'] = 'l'

        self.run_plotit(configuration)

        self.compare_images(
                os.path.join(self.output_folder.name, 'histo1.pdf'),
                get_golden_file('default_configuration_multi_stacks_line_type.pdf')
                )
