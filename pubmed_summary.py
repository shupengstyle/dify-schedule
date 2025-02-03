import requests
import smtplib
import os
import sys
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from bs4 import BeautifulSoup
from Bio import Entrez
from googletrans import Translator
from sumy.parsers.plaintext import PlaintextParser
from sumy.nlp.tokenizers import Tokenizer
from sumy.summarizers.lsa import LsaSummarizer

# 配置 PubMed API
Entrez.email = "your-email@example.com"  # 替换为你的邮箱

def summarize_text(text, sentences=3):
    """使用 LSA 进行文本摘要"""
    parser = PlaintextParser.from_string(text, Tokenizer("english"))
    summarizer = LsaSummarizer()
    summary = summarizer(parser.document, sentences)
    return "\n".join([str(sentence) for sentence in summary])

def fetch_pubmed_articles(query="COVID-19", max_results=5):
    """获取 PubMed 文章摘要或全文，并生成总结"""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="pub_date")
    record = Entrez.read(handle)
    handle.close()
    
    article_ids = record["IdList"]
    
    if not article_ids:
        return "No articles found."
    
    translator = Translator()
    articles = []
    
    for pmid in article_ids:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text")
        xml_data = handle.read()
        handle.close()
        
        soup = BeautifulSoup(xml_data, "lxml")
        abstract = soup.find("abstracttext")
        title = soup.find("articletitle").text if soup.find("articletitle") else "No title"
        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        
        # 尝试获取 PMC 免费全文链接
        pmc_link = None
        for article_id in soup.find_all("articleid"):
            if article_id.get("idtype") == "pmc":
                pmc_link = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{article_id.text}/"
                break

        summary_text = ""
        
        if pmc_link:
            try:
                full_text = requests.get(pmc_link).text
                soup = BeautifulSoup(full_text, "html.parser")
                sections = {s.text.lower(): s.find_next("p").text for s in soup.find_all("h2")}
                
                background = sections.get("background", "")
                methods = sections.get("methods", "")
                conclusion = sections.get("conclusion", "")
                
                summary_text = f"**背景:** {background}\n\n**方法:** {methods}\n\n**结论:** {conclusion}"
            except Exception as e:
                summary_text = f"无法解析全文: {e}"
        else:
            if abstract:
                summary_text = summarize_text(abstract.text)
            else:
                summary_text = "无摘要可用"

        translated_summary = translator.translate(summary_text, src="en", dest="zh-cn").text
        
        articles.append(f"PMID: {pmid}\n标题: {title}\n原文链接: {link}\n全文链接: {pmc_link if pmc_link else '无'}\n\n**摘要总结:**\n{translated_summary}\n{'-'*40}")
    
    with open("pubmed_articles.txt", "w", encoding="utf-8") as f:
        f.write("\n\n".join(articles))

    return "Articles summarized successfully."

def send_email():
    """发送邮件"""
    EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")
    EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")
    SMTP_SERVER = os.getenv("EMAIL_SMTP_SERVER")
    SMTP_PORT = int(os.getenv("EMAIL_SMTP_PORT", 587))

    with open("pubmed_articles.txt", "r", encoding="utf-8") as f:
        articles_content = f.read()

    subject = "最新 PubMed 文献 (含全文/摘要总结)"
    msg = MIMEMultipart()
    msg["From"] = EMAIL_ADDRESS
    msg["To"] = EMAIL_ADDRESS
    msg["Subject"] = subject
    msg.attach(MIMEText(articles_content, "plain"))

    server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
    server.starttls()
    server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)
    server.sendmail(EMAIL_ADDRESS, EMAIL_ADDRESS, msg.as_string())
    server.quit()

if __name__ == "__main__":
    if "--send" in sys.argv:
        send_email()
    else:
        print(fetch_pubmed_articles())
