# Deploying LocusCompare

The site is two things, both fronted by **Caddy** (automatic HTTPS):

- **locuscompare.com** → a static landing page (`../landing/`) linking to the four LocusCompare tools.
- **app.locuscompare.com** → the interactive single-locus Shiny app (`../`, built by `../Dockerfile`).

The database is **external** — the app connects to the shared LocusCompare2 reference DB
(`colotool` at `54.254.162.217:31987`), so there is no database container to run.

---

## 1. Prerequisites

- An **AWS EC2** instance running Ubuntu (a `t3.small`, 2 vCPU / 2 GB, is plenty — the app is
  light and the database is remote). Give it an **Elastic IP**.
- Security group: inbound **80** and **443** from anywhere, **22** from your admin IP only.
- **Docker Engine + Docker Compose plugin** installed.
- The `colotool` database firewall must **allow this EC2's public IP** to reach
  `54.254.162.217:31987` (confirm with the DB owner; it currently accepts our test host).

> There is already an older EC2 (`ec2-35-85-198-214.us-west-2`) that locuscompare.com forwards
> to and whose Shiny app returns a 500. You can either redeploy onto that instance or provision a
> fresh one; either way the steps below are the same.

## 2. Deploy

```bash
git clone https://github.com/boxiangliu/locuscompare.git
cd locuscompare
git checkout single-locus-redeploy          # until merged to master

cp .env.example .env
# edit .env and set DB_PASSWORD (the read-only colotool password)

docker compose up -d --build
```

This builds the app image (installs R + the latest `locuscomparer` from GitHub), starts the
Shiny container, and starts Caddy. Caddy will obtain Let's Encrypt certificates automatically
**once DNS points at this host and ports 80/443 are reachable** (next step).

Check status: `docker compose ps` and `docker compose logs -f caddy app`.

## 3. DNS (Namecheap)

locuscompare.com is managed at Namecheap and currently uses a **URL-forward** to the old EC2.
Replace it with A records pointing at the new Elastic IP:

| Host   | Type | Value            |
|--------|------|------------------|
| `@`    | A    | `<ELASTIC_IP>`   |
| `www`  | A    | `<ELASTIC_IP>`   |
| `app`  | A    | `<ELASTIC_IP>`   |

Remove the existing "URL Redirect / Forwarding" record for `@`/`www`. DNS propagation can take a
few minutes; Caddy issues certificates on the first HTTPS request after the names resolve.

## 4. Verify

- `https://locuscompare.com` → landing page with four working cards, valid certificate.
- `https://app.locuscompare.com` → the app. Click **Go to website → Single Locus**, download the
  two **Example file** links, upload them (labels e.g. "eQTL" / "GWAS"), enter reference SNP
  `rs9349379`, click **Plot LocusCompare** → the LocusCompare scatter + two LocusZoom panels
  render, colored by LD; switching **Population** recolors; **Download** returns a results ZIP.

## 5. Operations & notes

- **Update the app:** `git pull && docker compose up -d --build`.
- **Logs:** `docker compose logs -f app` (Shiny) / `caddy` (TLS/proxy).
- **Secrets:** `.env` is gitignored — never commit it. Rotate the old exposed credentials
  (`config/locuscompare-gmail-api-key.json`, `maintenance/*.pem`, the old Gmail app password);
  they are unused by this app and are already gitignored, but should be revoked and removed.
- **Database reliability:** the app depends on the external `colotool` host being reachable. If
  it becomes unavailable, plots that need position/LD lookups will fail. If you later want a
  self-contained deployment, dump the reference tables (`tkg_p3v5a_*`, `gencode_v19_gtex_v6p_hg38`)
  into a DB you control and repoint `.env` / the `locuscomparer` config.
- **hg19 gene lookup:** the read-only user cannot read `gencode_v19_gtex_v6p` (hg19), so
  gene-name region selection works on **hg38** only; SNP and coordinate selection work on both.
  Ask the DB owner for `GRANT SELECT ON colotool.gencode_v19_gtex_v6p` to enable hg19 gene lookup.

## 6. Local testing (no AWS)

```bash
cp .env.example .env   # set DB_PASSWORD
docker compose up -d --build app
# app is on the container's 3838; add a temporary "ports: 3838:3838" to the app service to hit it
```

Or run directly with R (needs the packages in the Dockerfile installed locally):

```bash
DB_PASSWORD=... APP_HOME=$PWD R -e 'shiny::runApp(".", port=8788, host="127.0.0.1")'
```

## 7. Optional: Terraform

If you prefer infrastructure-as-code, a minimal config would declare: `aws_instance`
(t3.small, Ubuntu 22.04 AMI), an `aws_eip` + association, an `aws_security_group`
(22/admin, 80+443/world), and `user_data` that installs Docker + Compose, clones this repo,
writes `.env`, and runs `docker compose up -d --build`. `terraform destroy` then gives a clean
teardown. This is optional — the manual steps above are sufficient.
