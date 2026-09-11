"""Read-only browser checks against an already running local studio."""
import json
from pathlib import Path
from playwright.sync_api import sync_playwright

with sync_playwright() as p:
    browser = p.chromium.launch(channel="chrome", headless=True, args=["--no-sandbox"])
    page = browser.new_page(viewport={"width": 1440, "height": 1000})
    reports = []
    for route in ["/", "/simulation", "/simulation/run", "/experiments", "/analysis"]:
        errors = []
        def on_error(error): errors.append(str(error))
        page.on("pageerror", on_error)
        response = page.goto("http://127.0.0.1:8092" + route, wait_until="networkidle")
        page.wait_for_timeout(1000)
        report = {"route": route, "status": response.status, "errors": errors,
                  "css": page.evaluate("document.styleSheets.length"),
                  "horizontal_overflow": page.evaluate("document.documentElement.scrollWidth > innerWidth + 2"),
                  "canvas_count": page.locator("canvas").count()}
        name = route.strip('/').replace('/', '-') or 'modeling'
        page.screenshot(path=f"/tmp/nanoplastic-{name}.png", full_page=True)
        reports.append(report)
        page.remove_listener("pageerror", on_error)
    browser.close()
    print(json.dumps(reports, indent=2))
    assert all(r['status'] == 200 and r['css'] and not r['errors'] and not r['horizontal_overflow'] for r in reports), reports
