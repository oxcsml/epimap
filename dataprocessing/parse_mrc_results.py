import json
from lxml import html
import pandas as pd
import requests


if __name__ == '__main__':
    page = requests.get('https://raw.githubusercontent.com/joshuablake/public-RTM-reports/master/2020-11-11.html')
    tree = html.fromstring(page.content)

    # R_t is under the following <div>:
    # r_t_div = tree.xpath('//div[@id="r_t"]')
    # and inside it, the data is stored as a json string as part of a widget.
    # The widget's identifier probably changes between weekly versions, and a
    # cleaner way is need to get the data (from underneath r_t_div):

    data = tree.xpath('//script[@data-for="htmlwidget-6b8dca356c43d597efec"]/text()')
    json_data = json.loads(data[0])

    annotations = [a['text'] for a in json_data['x']['layout']['annotations']]
    intervals = ['Median', 'Lower 95% CrI', 'Upper 95% CrI']

    df = pd.DataFrame(columns=['NHS Region', 'Date'] + intervals) \
           .set_index(['NHS Region', 'Date'])

    for i, element in enumerate(json_data['x']['data']):
        area = annotations[i // 3]
        dates = element['x']
        plot_annotations = element['hovertemplate']
        interval = intervals[0]
        if any(intervals[1] in s for s in plot_annotations):
            interval = intervals[1]
        if any(intervals[2] in s for s in plot_annotations):
            interval = intervals[2]

        # The json data in the html has a bug, which doesn't surface when
        # plotted. The 95% upper confidence interval contains the lower
        # confidence interval too in its time series (which is hidden in the
        # plot.
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
            df = df.append(df2)

    df = df.sort_index()
    df.to_csv('../data/comparisons/MRC-BSU-Cambridge-nowcasting-and-forecasting.csv')
