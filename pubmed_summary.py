import os
import requests
from Bio import Entrez
from bs4 import BeautifulSoup
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from googletrans import Translator

# 获取 PubMed API Key & Email 配置
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY")
EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")
SMTP_SERVER = os.getenv("EMAIL_SMTP_SERVER", "smtp.yeah.net")
SMTP_PORT = int(os.getenv("EMAIL_SMTP_PORT", 465))

# 配置 Entrez API
Entrez.email = EMAIL_ADDRESS
Entrez.api_key = PUBMED_API_KEY
translator = Translator()

def fetch_pubmed_articles(query="COVID-19", max_results=5):
    """获取 PubMed 文献，读取全文（如果可用），否则使用摘要"""
    
    # PubMed 搜索
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="pub_date", api_key=PUBMED_API_KEY)
    record = Entrez.read(handle)
    handle.close()

    article_ids = record["IdList"]
    if not article_ids:
        return "❌ 未找到相关文献。"

    articles = []

    for idx, pmid in enumerate(article_ids, start=1):
        # 获取文献详情
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text", api_key=PUBMED_API_KEY)
        xml_data = handle.read()
        handle.close()

        # 解析 XML
        soup = BeautifulSoup(xml_data, "lxml-xml")
        title = soup.find("articletitle").text if soup.find("articletitle") else "无标题"
        
        # 提取作者信息
        authors = soup.find_all("author")
        author_names = ', '.join([f"{author.find('lastname').text} {author.find('initials').text}" for author in authors if author.find('lastname')])
        
        # 提取期刊名称，使用大写的 <Journal> 标签
        journal_name = soup.find("Journal").find("isoabbreviation").text if soup.find("Journal").find("isoabbreviation") else "无杂志名称"
        
        # 获取发表年份
        pub_date = soup.find("pubdate")
        pub_year = pub_date.find("year").text if pub_date and pub_date.find("year") else "未知年份"
        
        # 摘要和 DOI
        abstract = soup.find("abstracttext")
        doi = soup.find("elocationid").text if soup.find("elocationid") else "无 DOI"
        full_text_url = f"https://doi.org/{doi}" if doi != "无 DOI" else None
        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        
        # 获取中文翻译
        translated_title = translator.translate(title, src="en", dest="zh-cn").text
        translated_abstract = translator.translate(abstract.text if abstract else "No abstract", src="en", dest="zh-cn").text

        # 格式化邮件内容
        article_content = f"""
{idx}. {title}. by {author_names} ({pub_year}) {journal_name}
 {translated_title}
PMID: {pmid} doi: {doi}
        """
        articles.append(article_content)

    # 邮件内容的文本格式
    text_content = "\n\n".join(articles)

    with open("pubmed_articles.txt", "w", encoding="utf-8") as f:
        f.write(text_content)

    return "✅ 文献摘要已获取并保存。"

def send_email():
    """发送邮件"""
    
    with open("pubmed_articles.txt", "r", encoding="utf-8") as f:
        articles_content = f.read()

    subject = "最新 PubMed 文献 (含全文/摘要总结)"
    msg = MIMEMultipart()
    msg["From"] = EMAIL_ADDRESS
    msg["To"] = EMAIL_ADDRESS
    msg["Subject"] = subject
    msg.attach(MIMEText(articles_content, "plain", "utf-8"))

    try:
        # 使用 SSL 连接 Yeah.net 邮件服务器
        server = smtplib.SMTP_SSL(SMTP_SERVER, SMTP_PORT)
        server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
        server.sendmail(EMAIL_ADDRESS, EMAIL_ADDRESS, msg.as_string())
        server.quit()
        print("📧 邮件发送成功！")
    except Exception as e:
        print(f"⚠️ 邮件发送失败: {e}")

if __name__ == "__main__":
    print(fetch_pubmed_articles())
    send_email()
