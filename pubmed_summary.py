import os
import requests
from Bio import Entrez
from bs4 import BeautifulSoup
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

# 获取 PubMed API Key & Email 配置
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY")
EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")
SMTP_SERVER = os.getenv("EMAIL_SMTP_SERVER", "smtp.yeah.net")
SMTP_PORT = int(os.getenv("EMAIL_SMTP_PORT", 465))

# 配置 Entrez API
Entrez.email = EMAIL_ADDRESS
Entrez.api_key = PUBMED_API_KEY

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

    for pmid in article_ids:
        # 获取文献详情
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="text", api_key=PUBMED_API_KEY)
        xml_data = handle.read()
        handle.close()

        # 解析 XML
        soup = BeautifulSoup(xml_data, "lxml")
        title = soup.find("articletitle").text if soup.find("articletitle") else "无标题"
        abstract = soup.find("abstracttext")
        full_text_url = None

        # 检查是否有全文链接
        for link in soup.find_all("elocationid"):
            if "doi" in link.text.lower():
                full_text_url = f"https://doi.org/{link.text}"

        link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

        if full_text_url:
            # 如果有全文链接，尝试获取全文内容
            summary_text = f"全文链接: {full_text_url}\n"
        else:
            # 否则使用摘要
            summary_text = abstract.text if abstract else "❌ 无摘要可用"

        articles.append(f"🔹 PMID: {pmid}\n📖 标题: {title}\n🔗 PubMed 链接: {link}\n📜 总结:\n{summary_text}\n{'-'*40}")

    with open("pubmed_articles.txt", "w", encoding="utf-8") as f:
        f.write("\n\n".join(articles))

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

