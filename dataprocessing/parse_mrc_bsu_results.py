
from html.parser import HTMLParser


import json
from lxml import html
import requests
import pandas as pd

class MyHTMLParser(HTMLParser):
    def handle_starttag(self, tag, attrs):
        print("Encountered a start tag:", tag)

    def handle_endtag(self, tag):
        print("Encountered an end tag :", tag)

    def handle_data(self, data):
        print("Encountered some data  :", data)


def p():
    page = requests.get('http://econpy.pythonanywhere.com/ex/001.html')
    with open('/Users/ulrich/Downloads/2020-11-11.html', 'r') as f:
        page = f.read()

    tree = html.fromstring(page)  #.content)

    r_t_div = tree.xpath('//div[@id="r_t"]')

    data = tree.xpath('//script[@data-for="htmlwidget-6b8dca356c43d597efec"]/text()')
    json_data = json.loads(data[0])

    # now write output to a file
    #twitterDataFile = open('/Users/ulrich/Downloads/2020-11-11.json', 'w')
    #twitterDataFile.write(json.dumps(j, indent=4, sort_keys=True))
    #twitterDataFile.close()

    annotations = [a['text'] for a in json_data['x']['layout']['annotations']]
    intervals = ['Median', 'Lower 95% CrI', 'Upper 95% CrI']

    df = pd.DataFrame(columns=['NHS Region', 'Date'] + intervals).set_index(['NHS Region', 'Date'])

    print(df)
    print(annotations)

    for i, element in enumerate(json_data['x']['data']):
        area = annotations[i // 3]
        dates = element['x']
        plot_annotations = element['hovertemplate']
        interval = intervals[0]
        if any(intervals[1] in s for s in plot_annotations):
            interval = intervals[1]
        if any(intervals[2] in s for s in plot_annotations):
            interval = intervals[2]

        print(interval)

        if interval == 'Upper 95% CrI':
            length = len(dates) // 2
        else:
            length = len(dates)

        data = {
            'NHS Region': [area for _ in range(length)],
            'Date': dates[-length:],
            interval: element['y'][-length:]
        }

        if i % 3 == 0:
            df2 = pd.DataFrame(data=data).set_index(['NHS Region', 'Date'])
        else:
            df2 = df2.join(pd.DataFrame(data=data).set_index(['NHS Region', 'Date']))

        if i % 3 == 2:
            # print(df2)
            df = df.append(df2)

    print(df)
    df.to_csv('/Users/ulrich/Downloads/MRC-BSU-Cambridge-nowcasting-and-forecasting.csv')
    #print(r_t_div)
    #<script type="application/json" data-for="htmlwidget-6b8dca356c43d597efec">
    #subtree = r_t_div.xpath('//script[]/text()')
    #print(subtree)

if __name__ == '__main__':
    p()
    # parser = MyHTMLParser()
    #parser.feed('<html><head><title>Test</title></head>'
    #        '<body><h1>Parse me!</h1></body></html>')

