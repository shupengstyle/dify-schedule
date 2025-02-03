import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from Bio import Entrez
from bs4 import BeautifulSoup
import google.generativeai as genai

# 环境变量配置
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY")
EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")
SMTP_SERVER = os.getenv("EMAIL_SMTP_SERVER", "smtp.yeah.net")
SMTP_PORT = int(os.getenv("EMAIL_SMTP_PORT", 465))
SEARCH_QUERY = os.getenv("SEARCH_QUERY")
MAX_RESULTS = int(os.getenv("MAX_RESULTS", 5))
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY") # 新增 Gemini API KEY

# 初始化配置
Entrez.email = EMAIL_ADDRESS
Entrez.api_key = PUBMED_API_KEY

# 配置 Gemini API
genai.configure(api_key=GOOGLE_API_KEY)
model = genai.GenerativeModel('gemini-pro')


def translate_text(text, target_language="zh-CN"):
    """使用 Gemini API 翻译文本"""
    if not text:
        return ""
    try:
        response = model.generate_content(f"Translate the following text to {target_language}: {text}")
        return response.text.strip() if response.text else text
    except Exception as e:
        print(f"Gemini API 翻译失败: {str(e)}")
        return text


def get_article_details(pmid):
    """获取单篇文献详细信息"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text")
        xml_data = handle.read()
        soup = BeautifulSoup(xml_data, "lxml-xml")

        # 解析元数据
        meta = {
            "title": soup.find("ArticleTitle").get_text() if soup.find("ArticleTitle") else "无标题",
            "journal": soup.find("ISOAbbreviation").get_text() if soup.find("ISOAbbreviation") else "未知杂志",
            "year": soup.find("Year").get_text() if soup.find("Year") else "未知年份",
            "doi": soup.find("ELocationID", {"EIdType": "doi"}).get_text() if soup.find("ELocationID", {"EIdType": "doi"}) else "无DOI"
        }

        # 处理作者信息
        authors = []
        for author in soup.find_all("Author"):
            lastname = author.find("LastName").get_text() if author.find("LastName") else ""
            initials = author.find("Initials").get_text() if author.find("Initials") else ""
            if lastname and initials:
                authors.append(f"{lastname} {initials}")

        # 处理摘要
        abstract = ""
        if abstract_section := soup.find("Abstract"):
            abstract = "\n".join([text.get_text() for text in abstract_section.find_all("AbstractText")])

        # 翻译处理
        translated_title = translate_text(meta["title"])
        translated_abstract = translate_text(abstract or "无摘要")

        return {
            "pmid": pmid,
            **meta,
            "authors": authors,
            "abstract": abstract,
            "translated_title": translated_title,
            "translated_abstract": translated_abstract,
            "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        }

    except Exception as e:
        print(f"获取文献 {pmid} 详情失败: {str(e)}")
        return None


def fetch_articles():
    """获取文献列表"""
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=SEARCH_QUERY,
            retmax=MAX_RESULTS,
            sort="pub_date",
            usehistory="y"
        )
        result = Entrez.read(handle)
        return result["IdList"]
    except Exception as e:
        print(f"文献搜索失败: {str(e)}")
        return []


def send_email(articles):
    """发送文献汇总邮件"""
    try:
        # 构建邮件内容
        msg = MIMEMultipart()
        msg["From"] = EMAIL_ADDRESS
        msg["To"] = EMAIL_ADDRESS
        msg["Subject"] = f"PubMed文献更新 - {SEARCH_QUERY}"

        # HTML内容模板
        html_content = f"""
        <html>
            <body>
                <h2>最新 {len(articles)} 篇文献</h2>
                <p>搜索关键词: {SEARCH_QUERY}</p>
        """

        for idx, article in enumerate(articles, 1):
            html_content += f"""
            <div style="margin-bottom: 30px; border-bottom: 1px solid #eee;">
                <h3>{idx}. {article['title']}</h3>
                <p><b>中文标题:</b> {article['translated_title']}</p>
                <p><b>作者:</b> {', '.join(article['authors'])}</p>
                <p><b>期刊:</b> {article['journal']} ({article['year']})</p>
                <p><b>摘要:</b><br>{article['abstract'] or '无摘要'}</p>
                <p><b>中文摘要:</b><br>{article['translated_abstract']}</p>
                <p>
                    <a href="{article['link']}">PubMed链接</a> | 
                    <a href="https://doi.org/{article['doi']}">全文链接</a>
                </p>
            </div>
            """

        html_content += "</body></html>"
        msg.attach(MIMEText(html_content, "html", "utf-8"))

        # 发送邮件
        with smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT) as server:
            server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
            server.send_message(msg)
            print("✅ 邮件发送成功")

    except Exception as e:
        print(f"邮件发送失败: {str(e)}")


if __name__ == "__main__":
    print("🚀 开始获取文献...")
    article_ids = fetch_articles()

    if not article_ids:
        print("❌ 未找到相关文献")
        exit()

    articles = []
    for pmid in article_ids:
        if article := get_article_details(pmid):
            articles.append(article)

    if articles:
        send_email(articles)
    else:
        print("❌ 未获取到有效文献信息")
