import requests
from lxml import etree
import xml.etree.ElementTree as ET
from collections import namedtuple

DBLP_BASE_URL = 'http://dblp.org/'
DBLP_AUTHOR_SEARCH_URL = DBLP_BASE_URL + 'search/author'

DBLP_PERSON_URL = DBLP_BASE_URL + 'pers/xk/{urlpt}'
DBLP_PUBLICATION_URL = DBLP_BASE_URL + 'rec/bibtex/{key}.xml'
DBLP_ALLPUB_URL = DBLP_BASE_URL+'pid/{pid}.xml'

class LazyAPIData(object):
    def __init__(self, lazy_attrs):
        self.lazy_attrs = set(lazy_attrs)
        self.data = None

    def __getattr__(self, key):
        if key in self.lazy_attrs:
            if self.data is None:
                self.load_data()
            return self.data[key]
        raise AttributeError(key)

    def load_data(self):
        pass

class Author(LazyAPIData):
    """
    Represents a DBLP author. All data but the author's key is lazily loaded.
    Fields that aren't provided by the underlying XML are None.

    Attributes:
    name - the author's primary name record
    publications - a list of lazy-loaded Publications results by this author
    homepages - a list of author homepage URLs
    homonyms - a list of author aliases
    """
    def __init__(self, urlpt):
        self.urlpt = urlpt
        self.xml = None
        super(Author, self).__init__(['name','pid','publications'])

    def load_data(self):
        resp = requests.get(DBLP_PERSON_URL.format(urlpt=self.urlpt))
        # TODO error handling
        xml = resp.content
        self.xml = xml
        
        root = etree.fromstring(xml)
        # print(DBLP_PUBLICATION_URL.format(key=self.key))

        pid = root.attrib['pid']
        resp_pub = requests.get(DBLP_ALLPUB_URL.format(pid=pid))
        self.publication_xml = resp_pub.content  # Store the fetched publication XML
        
        # XML 문자열을 파싱합니다.
        root = ET.fromstring(self.publication_xml)

        # 모든 서지 정보를 담을 리스트 생성
        publications = []

        # 모든 논문(r)을 반복하며 정보를 추출하여 리스트에 추가합니다.
        for r in root.findall('.//*'):
            if r.tag in ['article', 'inproceedings', 'proceedings', 'incollection', 'book', 'phdthesis', 'mastersthesis','data']:
                publication = {}
                
                publication['Author'] = root.attrib['name']
                publication['year'] = int(r.find('year').text)
                journal = r.find('journal')
                booktitle = r.find('booktitle')
                

                if journal is not None:
                    venue = journal.text
                elif booktitle is not None:
                    venue = booktitle.text
                else:
                    venue = "N/A"
                
                publication['venue'] = venue
                publication['title'] = r.find('title').text
                publications.append(publication)

        data = {
            'name':root.attrib['name'],
            'pid':pid,
            'publications':publications
        }

        self.data = data


def first_or_none(seq):
    try:
        return next(iter(seq))
    except StopIteration:
        pass

Publisher = namedtuple('Publisher', ['name', 'href'])
Series = namedtuple('Series', ['text','href'])
Citation = namedtuple('Citation', ['reference','label'])

class Publication(LazyAPIData):
    """
    Represents a DBLP publication- eg, article, inproceedings, etc. All data but
    the key is lazily loaded. Fields that aren't provided by the underlying XML
    are None.

    Attributes:
    type - the publication type, eg "article", "inproceedings", "proceedings",
    "incollection", "book", "phdthesis", "mastersthessis"
    sub_type - further type information, if provided- eg, "encyclopedia entry",
    "informal publication", "survey"
    title - the title of the work
    authors - a list of author names
    journal - the journal the work was published in, if applicable
    volume - the volume, if applicable
    number - the number, if applicable
    chapter - the chapter, if this work is part of a book or otherwise
    applicable
    pages - the page numbers of the work, if applicable
    isbn - the ISBN for works that have them
    ee - an ee URL
    crossref - a crossrel relative URL
    publisher - the publisher, returned as a (name, href) named tuple
    citations - a list of (text, label) named tuples representing cited works
    series - a (text, href) named tuple describing the containing series, if
    applicable
    """
    def __init__(self, key):
        self.key = key
        self.xml = None
        super(Publication, self).__init__( ['type', 'sub_type', 'mdate',
                'authors', 'editors', 'title', 'year', 'month', 'journal',
                'volume', 'number', 'chapter', 'pages', 'ee', 'isbn', 'url',
                'booktitle', 'crossref', 'publisher', 'school', 'citations',
                'series'])

    def load_data(self):
        resp = requests.get(DBLP_PUBLICATION_URL.format(key=self.key))
        
        xml = resp.content
        self.xml = xml
        root = etree.fromstring(xml)
        publication = first_or_none(root.xpath('/dblp/*[1]'))
        if publication is None:
            raise ValueError
        data = {
            # 'type':publication.tag,
            # 'sub_type':publication.attrib.get('publtype', None),
            # 'mdate':publication.attrib.get('mdate', None),
            # 'authors':publication.xpath('author/text()'),
            # 'editors':publication.xpath('editor/text()'),
            'title':first_or_none(publication.xpath('title/text()')),
            'year':int(first_or_none(publication.xpath('year/text()'))),
            # 'month':first_or_none(publication.xpath('month/text()')),
            'journal':first_or_none(publication.xpath('journal/text()')),
            # 'volume':first_or_none(publication.xpath('volume/text()')),
            # 'number':first_or_none(publication.xpath('number/text()')),
            # 'chapter':first_or_none(publication.xpath('chapter/text()')),
            # 'pages':first_or_none(publication.xpath('pages/text()')),
            # 'ee':first_or_none(publication.xpath('ee/text()')),
            # 'isbn':first_or_none(publication.xpath('isbn/text()')),
            # 'url':first_or_none(publication.xpath('url/text()')),
            'booktitle':first_or_none(publication.xpath('booktitle/text()'))#,
            # 'crossref':first_or_none(publication.xpath('crossref/text()')),
            # 'publisher':first_or_none(publication.xpath('publisher/text()')),
            # 'school':first_or_none(publication.xpath('school/text()')),
            # 'citations':[Citation(c.text, c.attrib.get('label',None))
            #              for c in publication.xpath('cite') if c.text != '...'],
            # 'series':first_or_none(Series(s.text, s.attrib.get('href', None))
            #           for s in publication.xpath('series'))
        }

        self.data = data

def search(author_str):
    resp = requests.get(DBLP_AUTHOR_SEARCH_URL, params={'xauthor':author_str})
    #TODO error handling
    root = etree.fromstring(resp.content)
    return [Author(urlpt) for urlpt in root.xpath('/authors/author/@urlpt')]
